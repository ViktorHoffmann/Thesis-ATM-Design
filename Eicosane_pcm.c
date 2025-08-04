//Modified UDF of the original source: https://akamcae.com/tutorials/phase-change-material-simulation-in-ansys-fluent/
#include "udf.h"
#include "mem.h"
#include "acc_profile_array.c"

//n-eicosane constant properties in solid phase
#define Ros_pcm 910.0 //Nazarychev-2022
#define Cps_pcm 2132.4 //NIST
#define Ks_pcm 0.4248 //Stryker-1990

//n-eicosane constant properties in fluid phase
#define Rol_pcm 769.0 //Nazarychev-2022
#define Cpl_pcm 2350.05 //NIST
#define Kl_pcm 0.1505 //Benbrika-2020

//thermal expansion coefficient
#define TEC 0.0009 //Benbrika-2020

//solidus and liquidus temperatures of n-eicosane
#define Ts 309.0 //NIST
#define Tl 311.0 //NIST

//reference temperature for Boussinesq's approximation
#define Tr 310.0
//please set the Tref in Fluent like Tr

//density of PCM
DEFINE_PROPERTY(Ro_var_PCM,cell,thread)
{
	double Gama, Ro_pcm;
	#if !RP_HOST
		Gama=C_LIQF(cell,thread);
		Ro_pcm=(1-Gama)*Ros_pcm+Gama*Rol_pcm;
	#endif
	return Ro_pcm;
}

DEFINE_SPECIFIC_HEAT(Cp_var_PCM,T,Tref,h,yi)
{
	double Gama, Cp_pcm;
	#if !RP_HOST
		if (T<Ts) { Cp_pcm=Cps_pcm; } else if (T>=Ts&&T<=Tl)
		{
			Gama=(T-Ts)/(Tl-Ts);
			Cp_pcm=((1-Gama)*Ros_pcm*Cps_pcm+Gama*Rol_pcm*Cpl_pcm)/((1-Gama)*Ros_pcm+Gama*Rol_pcm);
		}
		else
		{
			Cp_pcm=Cpl_pcm;
		}
		*h=Cp_pcm*(T-Tref);
	#endif
	return Cp_pcm;
}

//thermal conductivity of PCM
DEFINE_PROPERTY(K_var_PCM,cell,thread)
{
	double Gama, K_pcm;
	#if !RP_HOST
		Gama=C_LIQF(cell,thread);
		K_pcm=(1-Gama)*Ks_pcm+Gama*Kl_pcm;
	#endif
	return K_pcm;
}

//dynamic viscosity of PCM 
DEFINE_PROPERTY(Mu_var_PCM,cell,thread)
{
	double Temp,Mu_pcm;
	#if !RP_HOST
		Temp=C_T(cell,thread);
		Mu_pcm=(9*pow(10.,-4)*pow(Temp,2)-0.6529*Temp+119.94)*pow(10.,-3);
	#endif
	return Mu_pcm;
}

//Z-momentum source
DEFINE_SOURCE(Boussinesq_momentum_source, cell, thread, dS, eqn)
{
    double Temp = C_T(cell, thread);
    double acc = 9.81;  // default fallback

    // Get current simulation time
    real current_time = CURRENT_TIME;

    // Look up acceleration using linear interpolation
    int i;
    for (i = 0; i < NUM_ACC_POINTS - 1; ++i) {
        if (current_time >= acc_profile[i][0] && current_time <= acc_profile[i + 1][0]) {
            double t1 = acc_profile[i][0];
            double t2 = acc_profile[i + 1][0];
            double a1 = acc_profile[i][1];
            double a2 = acc_profile[i + 1][1];

            // Linear interpolation
            acc = a1 + (a2 - a1) * (current_time - t1) / (t2 - t1);
            break;
        }
    }

    // If outside range, use last known value
    if (current_time > acc_profile[NUM_ACC_POINTS - 1][0]) {
        acc = acc_profile[NUM_ACC_POINTS - 1][1];
    }

    double source = -Rol_pcm * acc * TEC * (Temp - Tr);
    dS[eqn] = -Rol_pcm * acc * TEC;

    return source;
}
