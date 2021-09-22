#pragma once

#include <vector>
#include <stdlib.h>
#include <iostream>
#include "functions.h" 		//for Bmag_dipole in Particles::u_and_p
#include "struct_Species.h"	//only for mass of electron
#include "constants.h" 		//only for speed of light(and L_shell, Re ->bounce period)
#include "common.h"

//Declaration of struct.

//Struct for particle's state throughout the iterations.
struct Particles
{ 				
	
	//Member function to push back initial gyrophases.
	void set_eta(real initial_eta);
	
	//Member function to push back initial equatorial pitch angles.
	void set_aeq(real initial_aeq);

	//Member function to push back initial latitudes.
	void set_lamda(real initial_lamda);
	
	void calculations0(real Beq0, real initial_lamda,real initial_zeta, real initial_time, real initial_aeq, real initial_Ekev);

	//Member function to push_back new state.
	void update_state(real new_aeq, real new_alpha, real new_lamda, real new_time);

	//Member variables.
	std::vector<real> lamda , zeta, uper , upar, ppar, pper, alpha, aeq, eta, M_adiabatic, time; 
};  	