#pragma once

#include <vector>
#include <stdlib.h>
#include "functions.h" //for Bmag_dipole in Particles::u_and_p
#include "struct_Species.h"	//only for mass of electron
#include "constants.h" 		//only for speed of light
#include "common.h"

//Declaration of struct.

//Struct for particle's state throughout the iterations.
struct Particles
{ 				
	
	//Constructor.
	Particles();
	
	//Member function to push back initial gyrophases.
	void set_eta(real initial_eta);
	
	//Member function to push back initial equatorial pitch angles.
	void set_aeq(real initial_aeq);
	
	void calculations0(real Beq0, real initial_lamda,real initial_zeta, real initial_time, real initial_M_adiabatic, real initial_aeq, real initial_Ekev);

	//Member function to push_back new state.
	void update_state(real new_lamda, real new_zeta, real new_uper, real new_upar, real new_ppar, real new_pper, real new_alpha, real new_alpha2, real new_aeq, real new_aeq2, real new_eta, real new_M_adiabatic, real new_time);

	//Member variables.
	std::vector<real> lamda , zeta, uper , upar, ppar, pper, alpha, alpha2, aeq, aeq2, eta, M_adiabatic, time; 
};  	