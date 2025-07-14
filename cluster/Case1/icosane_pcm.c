#include "udf.h"
#include "mem.h"

//n-eicosane constant properties in solid phase
#define Ros_pcm 910.0
#define Cps_pcm 1926.0
#define Ks_pcm 0.423

//n-eicosane constant properties in fluid phase
#define Rol_pcm 769.0
#define Cpl_pcm 2400.0
#define Kl_pcm 0.146

//thermal expansion coefficient
#define TEC 0.0008161

//solidus and liquidus temperatures of n-eicosane
#define Ts 309.15
#define Tl 311.15

//reference temperature for Boussinesq's approximation
#define Tr 310.15 //please select based on your problem

//please set the Tref in Fluent like Tr

//density of PCM
DEFINE_PROPERTY(Ro_var_PCM,cell,thread)
{
	double Gama, Ro_pcm;
	Gama=C_LIQF(cell,thread);
	Ro_pcm=(1-Gama)*Ros_pcm+Gama*Rol_pcm;
	return Ro_pcm;
}

DEFINE_SPECIFIC_HEAT(Cp_var_PCM,T,Tref,h,yi)
{
	double Gama, Cp_pcm;
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
	return Cp_pcm;
}

//thermal conductivity of PCM
DEFINE_PROPERTY(K_var_PCM,cell,thread)
{
	double Gama, K_pcm;
	Gama=C_LIQF(cell,thread);
	K_pcm=(1-Gama)*Ks_pcm+Gama*Kl_pcm;
	return K_pcm;
}

//dynamic viscosity of PCM 
DEFINE_PROPERTY(Mu_var_PCM,cell,thread)
{
	double Temp,Mu_pcm;
	Temp=C_T(cell,thread);
	Mu_pcm=(9*pow(10.,-4)*pow(Temp,2)-0.6529*Temp+119.94)*pow(10.,-3);
	return Mu_pcm;
}

//Z-momentum source
DEFINE_SOURCE(Boussinesq_momentum_source,cell,thread,dS,eqn)
{
	double Temp, source;
	Temp=C_T(cell,thread);
	source=-Rol_pcm*9.81*TEC*(Temp-Tr); //negative for -Z down
	dS[eqn]=-Rol_pcm*9.81*TEC; //negative for -Z down
	return source;
}