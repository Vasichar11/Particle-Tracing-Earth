#pragma once
#include <filesystem>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
#include "common.h"
#include "constants.h"
#include "read_hdf5.h"
#include "struct_Particles.h"
//HDF5 I/O
#include <highfive/H5File.hpp>          
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>


struct SetupArgs{
    real t_nowpi;
    real t_wpi;
    real t; // Total simulation time
    int64_t Nsteps_wpi;
    int64_t Nsteps_nowpi;
    bool use_bell_equations;
    bool use_li_equations;
    bool wpi;
    bool nowpi;
};

void InputArguments(int argc, char** argv, SetupArgs& setupArgs);
std::filesystem::path SelectFile(const std::filesystem::path directory, const std::string& extension, const std::string& prefix);
void readDstr(const std::string& dstrFilepath, std::vector<Particles>& dstr);
