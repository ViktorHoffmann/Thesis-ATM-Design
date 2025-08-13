//Modified UDF of the original source: https://akamcae.com/tutorials/phase-change-material-simulation-in-ansys-fluent/
#include "udf.h"
#include "mem.h"

//n-eicosane constant properties in solid phase
#define Ros_pcm 910.0
#define Cps_pcm 2132.4
#define Ks_pcm 0.4248

//n-eicosane constant properties in fluid phase
#define Rol_pcm 769.0
#define Cpl_pcm 2350.05
#define Kl_pcm 0.1505

//thermal expansion coefficient
#define TEC 0.0009

//solidus and liquidus temperatures of n-eicosane
#define Ts 309.0
#define Tl 311.0

//reference temperature for Boussinesq's approximation
#define Tr 310.0		//Fluent Tref must be equal to Tr

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

//dynamic viscosity of PCM with fit
DEFINE_PROPERTY(Mu_var_PCM,cell,thread)
{
	double Temp,Mu_pcm;
	#if !RP_HOST
		Temp=C_T(cell,thread);
		Mu_pcm=(9*pow(10.,-4)*pow(Temp,2)-0.6529*Temp+119.94)*pow(10.,-3);
	#endif
	return Mu_pcm;
}

/*...*/

//Y-momentum source
DEFINE_SOURCE(Boussinesq_momentum_source,cell,thread,dS,eqn)
{
	double Temp, source, acc;
	Temp=C_T(cell,thread);

	double t = CURRENT_TIME;

	if (t < 20)
		acc = 34.81;
	else if (t < 50)
		acc = 109.81;
	else if (t < 150)
		acc = 19.62;
	else
		acc = 9.81;

	source=-Rol_pcm*acc*TEC*(Temp-Tr);  //negative for -Y down
	dS[eqn]=-Rol_pcm*acc*TEC; 			//negative for -Y down
	return source;
}