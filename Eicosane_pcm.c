// Modified UDF with transient acceleration from profile
#include "udf.h"
#include "mem.h"
#include "profile.h"

// === Constants ===

// n-eicosane solid phase
#define Ros_pcm 910.0
#define Cps_pcm 2132.4
#define Ks_pcm 0.4248

// n-eicosane liquid phase
#define Rol_pcm 769.0
#define Cpl_pcm 2350.05
#define Kl_pcm 0.1505

#define TEC 0.0009         // Thermal expansion coefficient
#define Ts 309.0           // Solidus temperature
#define Tl 311.0           // Liquidus temperature
#define Tr 310.0           // Reference temperature

// === Density ===
DEFINE_PROPERTY(Ro_var_PCM, cell, thread)
{
    double Gama, Ro_pcm;
    #if !RP_HOST
        Gama = C_LIQF(cell, thread);
        Ro_pcm = (1 - Gama) * Ros_pcm + Gama * Rol_pcm;
    #endif
    return Ro_pcm;
}

// === Specific Heat ===
DEFINE_SPECIFIC_HEAT(Cp_var_PCM, T, Tref, h, yi)
{
    double Gama, Cp_pcm;
    #if !RP_HOST
        if (T < Ts) {
            Cp_pcm = Cps_pcm;
        } else if (T <= Tl) {
            Gama = (T - Ts) / (Tl - Ts);
            Cp_pcm = ((1 - Gama) * Ros_pcm * Cps_pcm + Gama * Rol_pcm * Cpl_pcm) /
                     ((1 - Gama) * Ros_pcm + Gama * Rol_pcm);
        } else {
            Cp_pcm = Cpl_pcm;
        }
        *h = Cp_pcm * (T - Tref);
    #endif
    return Cp_pcm;
}

// === Thermal Conductivity ===
DEFINE_PROPERTY(K_var_PCM, cell, thread)
{
    double Gama, K_pcm;
    #if !RP_HOST
        Gama = C_LIQF(cell, thread);
        K_pcm = (1 - Gama) * Ks_pcm + Gama * Kl_pcm;
    #endif
    return K_pcm;
}

// === Dynamic Viscosity ===
DEFINE_PROPERTY(Mu_var_PCM, cell, thread)
{
    double Temp, Mu_pcm;
    #if !RP_HOST
        Temp = C_T(cell, thread);
        Mu_pcm = (9e-4 * pow(Temp, 2) - 0.6529 * Temp + 119.94) * 1e-3;
    #endif
    return Mu_pcm;
}

// === Boussinesq Momentum Source with Time-Varying Acceleration ===
static Profile *acc_profile = NULL;

DEFINE_SOURCE(Boussinesq_momentum_source, cell, thread, dS, eqn)
{
    real acc = 9.81; // default fallback value
    real Temp = C_T(cell, thread);
    real source;

    if (!acc_profile)
        acc_profile = Get_Profile("acc", thread->id);  // "acc" must match your profile name

    if (acc_profile)
        Profile_Eval(acc_profile, CURRENT_TIME, &acc, 1);

    source = -Rol_pcm * acc * TEC * (Temp - Tr);
    dS[eqn] = -Rol_pcm * acc * TEC;

    return source;
}