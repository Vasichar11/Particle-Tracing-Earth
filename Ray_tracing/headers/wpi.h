#pragma once
#include <vector>
#include <algorithm>   //min_element, max_element for std vector. Is it slow?
#include "read_hdf5.h" //read hdf5 ray data
#include "common.h"
#include "constants.h"
#include "functions.h"
#include "is_in_packet.h"
#include "time_rates.h"
#include "struct_Telescope.h"


#include <highfive/H5File.hpp>                                      
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;

//Function for WPI. Involves RK4 for Nsteps. 

void wpi(real p, real lamda, real alpha, real aeq, real ppar, real pper, real upar, real uper, real zeta, real M_adiabatic, real eta, real time, Telescope &ODPT);

