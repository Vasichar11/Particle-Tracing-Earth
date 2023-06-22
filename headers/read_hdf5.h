#pragma once
#include <string>
#include <vector>
#include <iostream>
#include "common.h"
#include <H5Cpp.h>

using namespace H5;

 std::vector<real> read_hdf5(const std::string& dataset_name, const std::string& file_name);
