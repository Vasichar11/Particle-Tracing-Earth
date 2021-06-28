#include"headers/struct_Species.h"


//Returns specie's lamda dependent density(estimation)
real Species::density(real l) {real clat=cos(l); 	real ns=(Constants::ne_0)*pow(clat,-4)*n_factor;  return ns; }

//Returns specie's plasma frequency squared
real Species::wps(real ns) { return ns*pow(charge,2)/(mass*(Constants::eps0)); } 

//Returns specie's cyclotron frequency
real Species::wc(real B0mag) { return (charge * B0mag)/mass; } 				   
