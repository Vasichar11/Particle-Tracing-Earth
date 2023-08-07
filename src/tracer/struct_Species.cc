#include"struct_Species.h"

Species::Species(real type_mass, real type_charge, real type_n_factor)
{
	mass     = type_mass;
	charge   = type_charge;
	n_factor = type_n_factor;
}

//Returns specie's latitude dependent density(estimation)
real Species::density(real l) {real clat=cos(l); 	real ns=(Universal::ne_0)*pow(clat,-4)*n_factor;  return ns; } //most called function 4x4xNstepsxPopulation
															//Denton et al., 2002
//Returns specie's plasma frequency squared
real Species::wps(real ns) { return ns*pow(charge,2)/(mass*(Universal::eps0)); } 

//Returns specie's cyclotron frequency
real Species::wc(real Bmag) { return (charge * Bmag)/mass; } 				   
