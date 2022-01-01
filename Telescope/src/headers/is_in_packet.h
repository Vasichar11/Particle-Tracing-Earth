#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include "common.h"
#include "constants.h"

//Check if the particle is within the wave's latitude.
//Return wave's index where there's matching latitude(WPI). Will be used for wave's data.
int is_in_packet(real min_lat, real max_lat, real lamda_tmp, int i, std::vector<real> &wave_lat);
