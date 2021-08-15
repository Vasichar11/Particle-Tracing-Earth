#pragma once
#include "common.h"
#include "constants.h"
#include <cmath>


//Struct for mass & charge data of species
//   and for computing densities & wps,wc
struct Species
{     	  	      				  
	//Returns specie's lamda dependent density(estimation)
	real density(real l);

	//Returns specie's plasma frequency squared
	real wps(real ns);

	//Returns specie's cyclotron frequency
	real wc(real B0mag);

	//Declare Member Variables
	real mass;		   	  	 	  
	real charge;
	real n_factor;        	 	  //For now, n_factor: e->1 , H->0.94, He->0.0054, O->0.006

};