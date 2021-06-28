#pragma once

#include <vector>
#include <stdlib.h>
#include "common.h"

//Declaration of struct.

//Struct for particle's state throughout the iterations.
struct Particles
{ 				
	
	//Constructor for initial state, called once when object is created.
	Particles(real first_lamda, real first_zeta, real first_upar, real first_uper, real first_ppar, real first_pper, real first_alpha, real first_alpha2, real first_aeq, real first_aeq2,real first_M_adiabatic,real first_time);
	
	//Member function to push back initial gyrophases.
	void set_eta(real eta_distribution);
	
	//Member function to push_back new state.
	void update_state(real new_lamda, real new_zeta, real new_uper, real new_upar, real new_ppar, real new_pper, real new_alpha, real new_alpha2, real new_aeq, real new_aeq2, real new_eta, real new_M_adiabatic, real new_time);

	//Member variables.
	std::vector<real> lamda , zeta, uper , upar, ppar, pper, alpha, alpha2, aeq, aeq2, eta, M_adiabatic, time; 
};  	
