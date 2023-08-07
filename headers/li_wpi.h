#pragma once
#include <vector>
#include <algorithm> // min_element, max_element for std vector.
#include <omp.h>
#include "read_hdf5.h" //read hdf5 ray data
#include "common.h"
#include "constants.h"
#include "functions.h"
#include "RK_estimations.h"
#include "is_in_packet.h"
#include "time_rates.h"
#include "struct_Telescope.h"
#include "struct_Particles.h"
#include "struct_Ray.h"

// Function for WPI. Involves RK4 for Nsteps. 
void li_wpi(const int64_t Nstep_wpi, const int p, Particles &single, Telescope &ODPT, Ray &ray);


