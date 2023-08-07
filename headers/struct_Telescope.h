#pragma once
#include "common.h"
#include <iostream>
#include <vector>

struct Telescope
{
	// Satellite's position parameters 
	real latitude_deg;
	real L_shell;						 
	
	// Vectors to store detected particles
	std::vector<real> latitude, alpha, aeq, eta, time; 
	std::vector<int> id;

	// Constructor - Initializes position of satellite - resizes vectors
	Telescope(real lat, real L);
	
	void resize(size_t size = 0);

	// Returns true if particle crossed satellite
	bool crossing(real p1_latitude, real p2_latitude, real p_L_shell);	

	// Store detected particles in satellite's memory
	void store(int id, real latitude, real aeq, real alpha, real time);

};	