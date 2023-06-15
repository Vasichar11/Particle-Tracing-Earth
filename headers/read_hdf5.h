#pragma once
#include <highfive/H5File.hpp>                                      
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <string>
#include <vector>
#include "common.h"

namespace h5 = HighFive;

 std::vector<real> read_hdf5(const std::string& dataset_name, const std::string& file_name);
