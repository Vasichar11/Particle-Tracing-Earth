#pragma once
#include "common.h"
#include "functions.h"
#include "time_rates.h"
#include "constants.h"
#include "struct_Species.h"
#include "struct_Telescope.h"
#include <vector>
#include <iostream>
#include <cmath>
    
//Function for particle's motion. Involves RK4 for Nsteps. 

void motion(int64_t track_pop, int p, real lamda, real alpha, real aeq, real ppar, real pper, real upar, real uper, real zeta, real M_adiabatic, real eta, real time, Telescope &ODPT);

