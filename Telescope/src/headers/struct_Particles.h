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
	void initialize(real eta0, real aeq0, real lamda0, real Ekin0, real zeta0, real time0);
		
	//Member function to save particle states.
	void escaping_state(int id,real new_lamda, real new_aeq, real new_alpha, real new_time);
	void negative_state(int id);
	void high_state(int id);
	void nan_state(int id);
	//If vectors for save_state are needed;
	void save_state(int id, real new_alpha, real new_deta_dt, real new_time);
	

	void lamda_domain(real aeq0);


	//Member variables.
	bool escaped,trapped,negative,nan,high; //Characterize particle
	//Starting noWPI states
	real lamda0, zeta0, uper0, upar0, ppar0, pper0, alpha0, aeq0, eta0, M_adiabatic0, Ekin0, time0;
	//Ending noWPI and Starting WPI states
	real lamda00, zeta00, uper00, upar00, ppar00, pper00, alpha00, aeq00, eta00, M_adiabatic00, Ekin00, time00;
	//If lost
	real lamda_lost, alpha_lost, aeq_lost, time_lost, id_lost;  //, ppar_lost, pper_lost, zeta_lost, uper_lost, upar_lost,  eta_lost, M_adiabatic_lost, Ekin_lost;
	//If negative P.A
	real id_neg;
	//If higher than 180 P.A
	real id_high;
	//If nan
	real id_nan;


	//For any saved state
	int saved_id;
	real saved_deta_dt;
	real saved_alpha;
	real saved_time;	

	//If vectors for save_state are needed;
	//std::vector<real> id, lamda , zeta, uper , upar, ppar, pper, alpha, aeq, eta, M_adiabatic, deta_dt, Ekin, time; 
	std::vector<real> id, alpha, deta_dt, time; 

};  	
