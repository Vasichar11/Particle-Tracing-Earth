#pragma once
#include "common.h"
#include <iostream>
#include <vector>

struct Telescope
{

	Telescope(real lat, real L_parameter); //Constructor. Initialize position of satellite.
	
	bool crossing(real p1_lamda, real p2_lamda, real p_L_shell);	//Returns true if particle crossed satellite.

	//Function to push back detected particles.
	void store(int id, real lamda, real zeta, real uper, real upar, real ppar, real pper, real alpha, real alpha2, real aeq, real aeq2, real eta, real M_adiabatic, real time);

	//void update_look_dir();   					//When called, update look direction.
	
	//--Satellite's position parameters--//
	real L_shell;						 //L_shell of the satelite.
	real latitude;

	//----Detector's View parameters----//
	//Omni directional particle telescope. Longitudinal direction is not yet considered.
	//real view;							//Detector's Angular coverage(degrees).
	//real angle;							//Detector's Angle of axis(degrees).
	//real sector_range;					//Sector's 	 Angular coverage(degrees).
	//int  sector_num = view/sector_range;//Number of Sectors.
	//int  sector;						//Viewing sector identifier.

	//----Detector's Energy channels----//
	//real Ebin_max;					//Energies it can detect.
	//real Ebin_min;		
	
	//--Detector's Temporal resolution--//
	//real tRes;				 		    //Minimum time to distinguish two events.
	
	//Vectors to store detected particles.
	std::vector<real> lamda , zeta, uper , upar, ppar, pper, alpha, alpha2, aeq, aeq2, eta, M_adiabatic, time; 
	std::vector<int> id;
};	