#pragma once

#include "common.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>



//Function. Read filename. Return Vector of pairs <column name, column vector values>
std::vector<std::pair<std::string, std::vector<double>>> rcsv(std::string filename);
