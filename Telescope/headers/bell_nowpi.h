#pragma once
#include "common.h"
#include "functions.h"
#include "time_rates.h"
#include "constants.h"
#include "struct_Species.h"
#include "struct_Telescope.h"
#include "struct_Particles.h"
#include <vector>
#include <iostream>
#include <cmath>
    
//Function for particle's adiabatic motion. Involves RK4 for Nsteps. 
void adiabatic_motion(int p, Particles &single, Telescope &ODPT);


