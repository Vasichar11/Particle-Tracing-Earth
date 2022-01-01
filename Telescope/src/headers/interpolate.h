#pragma once

#include "common.h"
#include <vector>

bool getLeftNeighbor(std::vector<real> data , real target,size_t& lNeighbor);
std::vector <real> interpolate(std::vector <real> time_new, std::vector <real> timef, std::vector <real> vector_csv) ;

