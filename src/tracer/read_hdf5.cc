#include "read_hdf5.h"


//Function to read HDF5 datasets.
 std::vector<real> read_hdf5(const std::string& dataset_name, const std::string& file_name)
 { 
    std::vector<real> read_data;
    try
    {
        // Open the HDF5 file and the dataset
        H5::H5File file_in(file_name, H5F_ACC_RDONLY);
        H5::DataSet dataset = file_in.openDataSet(dataset_name);

        // Get the datatype, the dataspece and the number of elements of the dataset
        H5::DataType data_type = dataset.getDataType();
        H5::DataSpace data_space = dataset.getSpace();
        hsize_t num_elements = data_space.getSimpleExtentNpoints();

        // Resize the vector to accommodate the data
        read_data.resize(num_elements);
        // Read the data from the dataset
        dataset.read(&read_data[0], data_type);

        // Close the dataset and the file
        dataset.close();
        file_in.close();
    }
    catch (H5::Exception& e) 
    {
        // Handle any exceptions that can occur during HDF5 operations
        std::cerr << "HDF5 exception: " << e.getDetailMsg() << std::endl;
    }

    return read_data;
}