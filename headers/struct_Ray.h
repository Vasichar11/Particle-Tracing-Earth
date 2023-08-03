#pragma once
#include <vector>
#include <string>
#include "common.h"
#include "read_hdf5.h"

struct Ray
{
    std::vector<real> lat_int;
    std::vector<real> kx_ray;
    std::vector<real> kz_ray;
    std::vector<real> kappa_ray;
    std::vector<real> Bzw;
    std::vector<real> Ezw;
    std::vector<real> Bw_ray;
    std::vector<real> w1;
    std::vector<real> w2;
    std::vector<real> R1;
    std::vector<real> R2;

    void readRay(const std::string& rayFilepath);
};