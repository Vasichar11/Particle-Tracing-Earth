#include "read_hdf5.h"


//Function to read HDF5 datasets.
 std::vector<real> read_hdf5(const std::string& dataset_name, const std::string& file_name)
 { 
     std::vector<real> read_data;

    h5::File file_in(file_name, h5::File::ReadOnly);

    h5::DataSet dataset = file_in.getDataSet(dataset_name);

    dataset.read(read_data);

    return read_data;
}