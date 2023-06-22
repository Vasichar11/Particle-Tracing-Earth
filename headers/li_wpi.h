#pragma once
#include <vector>
#include <algorithm>   //min_element, max_element for std vector. Is it slow?
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

//Function for WPI. Involves RK4 for Nsteps. 

void li_wpi(const int64_t Nstep_wpi, const int p, std::vector <real> &lat_int, const std::vector <real> &kx_ray, const std::vector <real> &kz_ray, const std::vector <real> &kappa_ray, const std::vector <real> &Bzw, const std::vector <real> &Ezw, const std::vector <real> &Bw_ray, const std::vector <real> &w1, const std::vector <real> &w2, const std::vector <real> &R1, const std::vector <real> &R2, Particles &single, Telescope &ODPT);


