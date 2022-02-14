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
	void initialize(real eta0, real aeq0, real alpha0, real lamda0, real Ekev0, real Blam0, real zeta0, real time0);
		
	//Member function to save particle states.
	void save_state(int id,real new_lamda, real new_alpha, real new_aeq, real new_time);



	//Member variables.
	bool escaped;
	//Initials
	real lamda_init, zeta_init, uper_init, upar_init, ppar_init, pper_init, alpha_init, aeq_init, eta_init, M_adiabatic_init, Ekin_init, time_init;
	//When end simulation type
	real lamda_end, zeta_end, uper_end, upar_end, ppar_end, pper_end, alpha_end, aeq_end, eta_end, M_adiabatic_end, Ekin_end, time_end;
	//When lost
	real lamda_lost, alpha_lost, aeq_lost, time_lost, id_lost;  //, ppar_lost, pper_lost, zeta_lost, uper_lost, upar_lost,  eta_lost, M_adiabatic_lost, Ekin_lost;

	//If vectors are needed.
	//std::vector<real> id, lamda , zeta, uper , upar, ppar, pper, alpha, aeq, eta, M_adiabatic, deta_dt, Ekin, time; 

};  	
