#pragma once
#include "common.h"
#include "functions.h"
#include "time_rates.h"
#include "RK_estimations.h"
#include "constants.h"
#include "struct_Species.h"
#include "struct_Telescope.h"
#include "struct_Particles.h"
#include <vector>
#include <iostream>
#include <cmath>

//Function for wave=particle interaction. Involves RK4 for Nsteps. 
void bell_wpi(int p, Particles &single, Telescope &ODPT);