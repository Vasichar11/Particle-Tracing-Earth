#include "struct_Ray.h"

void Ray::readRay(const std::string& rayFilepath) {
    lat_int = read_hdf5("lat_int", rayFilepath);
    kx_ray = read_hdf5("kx_ray", rayFilepath);
    kz_ray = read_hdf5("kz_ray", rayFilepath);
    kappa_ray = read_hdf5("kappa_ray", rayFilepath);
    Bzw = read_hdf5("Bzw", rayFilepath);
    Ezw = read_hdf5("Ezw", rayFilepath);
    Bw_ray = read_hdf5("Bw_ray", rayFilepath);
    w1 = read_hdf5("w1", rayFilepath);
    w2 = read_hdf5("w2", rayFilepath);
    R1 = read_hdf5("R1", rayFilepath);
    R2 = read_hdf5("R2", rayFilepath);
}

