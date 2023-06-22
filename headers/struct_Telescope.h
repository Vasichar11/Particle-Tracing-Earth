#pragma once
#include "common.h"
#include <iostream>
#include <vector>

struct Telescope
{

	Telescope(real latitude, real L_shell); //Constructor. Initialize position of satellite.
	
	bool crossing(real p1_latitude, real p2_latitude, real p_L_shell);	//Returns true if particle crossed satellite.

	//Function to push back detected particles.
	void store(int id, real latitude, real aeq, real alpha, real time);

	//void update_look_dir();   					//When called, update look direction.
	
	//--Satellite's position parameters--//
	real L_shell;						 //L_shell of the satelite.
	real latitude;
	
	//Vectors to store detected particles.
	std::vector<real> latitude , uper , upar, alpha, aeq, eta, time; 
	std::vector<int> id;
};	