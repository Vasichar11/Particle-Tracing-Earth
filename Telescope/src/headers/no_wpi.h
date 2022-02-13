#pragma once
#include "common.h"
#include "functions.h"
#include "time_rates.h"
#include "constants.h"
#include "struct_Species.h"
#include "RK_estimations.h"
#include "struct_Telescope.h"
#include "struct_Particles.h"
#include <vector>
#include <iostream>
#include <cmath>
    
//Function for particle's adiabatic motion. Involves RK4 for Nsteps. 
void no_wpi(int p, Particles &single, Telescope &ODPT, Particles &particle_state);


