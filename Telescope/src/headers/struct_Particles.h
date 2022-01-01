#pragma once

#include <vector>
#include <stdlib.h>
#include <iostream>
#include "functions.h" 		//for Bmag_dipole in Particles::u_and_p
#include "constants.h" 		//only for speed of light(and L_shell, Re ->bounce period)
#include "common.h"

//Declaration of struct.

//Struct for particle's state throughout the iterations.
struct Particles
{ 				
	//Member function to initialize particle population.
	void initialize(real eta0, real aeq0, real alpha0, real lamda0, real Ekev0, real Blam0, real zeta0, real deta_dt0, real time0);
		
	//Member function to push_back new state if needed.
	void save_state(real new_lamda, real new_alpha, real new_aeq, real new_ppar, real new_pper, real new_time);

	//Member variables.
	std::vector<real> lamda , zeta, uper , upar, ppar, pper, alpha, aeq, eta, M_adiabatic, deta_dt, Ekin, time; 
};  	
