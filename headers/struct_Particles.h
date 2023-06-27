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
	// Characterize particle
	bool escaped, trapped, negative, nan, high; 
	// Starting noWPI states
	real latitude0, zeta0, uper0, upar0, ppar0, pper0, alpha0, aeq0, eta0, M_adiabatic0, Ekin0, time0;
	// Ending noWPI and Starting WPI states
	real latitude_end, zeta_end, uper_end, upar_end, ppar_end, pper_end, alpha_end, aeq_end, eta_end, M_adiabatic_end, Ekin_end, time_end;
	// If lost
	real latitude_lost, alpha_lost, aeq_lost, time_lost, id_lost;  //, ppar_lost, pper_lost, zeta_lost, uper_lost, upar_lost,  eta_lost, M_adiabatic_lost, Ekin_lost;
	// If negative P.A
	real id_neg;
	// If higher than 180 P.A
	real id_high;
	// If nan
	real id_nan;
	// For any saved state
	int saved_id;
	real saved_max_dEkin, saved_maxEkin_time, saved_max_dPA, saved_maxdPA_time, saved_mindetadt_time, saved_min_detadt;
	// If vectors for save_state are needed
	//std::vector<real> id, latitude , zeta, uper , upar, ppar, pper, alpha, aeq, eta, M_adiabatic, deta_dt, Ekin, time; 
	std::vector<real> id, latitude, deta_dt, Ekin; 

	// Member function to initialize particle population
	void initialize(real eta0, real aeq0, real latitude0, real Ekin0, real zeta0, real time0);
	// Member function to save particles that meet certain conditions
	void escaping_state(int id,real new_latitude, real new_aeq, real new_alpha, real new_time);
	void negative_state(int id);
	void high_state(int id);
	void nan_state(int id);
	// Member function to save other particle states
	void save_state(int p, real max_dEkin, real maxEkin_time, real  max_dPA, real maxdPA_time, real min_detadt, real min_detadt_time);
	

	void latitude_domain(real aeq0);




};  	
