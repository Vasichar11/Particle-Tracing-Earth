#include "parameters.h"
namespace h5 = HighFive;


void InputArguments(int argc, char** argv, SetupArgs& setupArgs) {

    // Validate number of arguments
    if (argc < 3  || argc > 4) {
        std::cerr << "Invalid number of arguments. Usage: program_name <noWPI_time> <WPI_time> [-bell]\n"
                << "\n- Argument variables are the simulation times for the NoWPI and WPI interaction in seconds\n"
                << "- First, provide the time(seconds) of NoWPI simulation and then the time(seconds) of the WPI simulation\n"
                << "- The default WPI simulation is the one using the Li Formulas (li_wpi.cc)\n"
                << "- For WPI using the Bell Formulas (bell_wpi.cc), run with the option -bell as the last command line argument\n";
        std::terminate();
    }

    std::string arg1 = argv[1];
    std::string arg2 = argv[2];

    // Validate arg1 and arg2 as numbers
    try {
        real number1 = std::stof(arg1); // Exception if arguments are not numbers
        real number2 = std::stof(arg2);
        (void)number1;  // Intentionally unused to avoid warnings
        (void)number2;  
    } catch (const std::invalid_argument&) {
        std::cerr << "Invalid argument. Argument 1 and argument 2 must be number of seconds .\n"
                << "\nArgument variables are the simulation times for the NoWPI and WPI interaction in seconds.\n"
                << "First, provide the time(seconds) of NoWPI simulation and then the time(seconds) of the WPI simulation.\n"
                << std::endl;
        std::terminate();
    }

    bool BooleanBell = false;
    std::string bell = "-bell";
    if (argc == 4) {
        std::string arg3 = argv[3];
        if (arg3 != bell) {
            std::cout << "Invalid argument. Argument 3 may only be \"-bell\".\n"
                    << "The default WPI simulation is the one using the Li Formulas (li_wpi.cc)\n"
                    << "For WPI using the Bell Formulas (bell_wpi.cc), run with the option -bell as the last command line argument.\n"
                    << std::endl;
            std::terminate();
        }
        else {
            BooleanBell = true;
        }
    }
    
    // Assign times from the command line arguments
    setupArgs.t_nowpi = std::stof(argv[1]); // No WPI time from command line 
    setupArgs.t_wpi = std::stof(argv[2]); // WPI time from command line 
    setupArgs.t = setupArgs.t_nowpi + setupArgs.t_wpi;  // Total simulation time
    setupArgs.Nsteps_wpi  = setupArgs.t_wpi/Simulation::h; // WPI step count
    setupArgs.Nsteps_nowpi = setupArgs.t_nowpi/Simulation::h; // noWPI step count
    setupArgs.use_bell_equations = BooleanBell; // To use Bell equations

    if (setupArgs.t_wpi > 0) setupArgs.wpi = true; else setupArgs.wpi = false;
    if (setupArgs.t_nowpi > 0) setupArgs.nowpi = true; else setupArgs.nowpi = false;
    if (setupArgs.use_bell_equations) setupArgs.use_li_equations = false; else setupArgs.use_li_equations = true;

    std::cout<<"Loss cone angle, for Lshell: "<<Distribution::L_shell<<" is " << Simulation::alpha_lc*Universal::R2D<< " degrees"<<std::endl;
    std::cout<<"NoWPI Simulation: "<<setupArgs.t_nowpi<< "s"<<std::endl;
    std::cout<<"WPI Simulation: "<<setupArgs.t_wpi<< "s"<<std::endl;
    
    // Success!
}

// Function to prompt the user for valid file selections
std::filesystem::path SelectFile(const std::filesystem::path directory, const std::string& extension, const std::string& prefix) {
    std::vector<std::string> files;
    std::string selectedFilename;
    std::filesystem::path selectedFilepath;

    // Iterate over the directory and store files that match the distribution hdf5 files
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        const std::string filename = entry.path().filename().string();
        const std::string stem = entry.path().stem().string();

        if (entry.path().extension() == extension && stem.substr(0, prefix.size()) == prefix) {
            files.push_back(filename);
        }
    }

    // Check if no files were found
    if (files.empty()) {
        throw std::runtime_error("No matching files found. Please make sure you have created them!");
    }

    int selection;
    bool validSelection = false;

    // Loop until valid file selection
    do {
        // Display
        for (size_t i = 0; i < files.size(); ++i) {
            std::cout << i + 1 << ". " << files[i] << '\n';
        }
        
        // Prompt selection
        std::cout << "Enter the number of the file you want to select: ";
        if (std::cin >> selection) {
            if (selection >= 1 && selection <= static_cast<int>(files.size())) {
                validSelection = true;
                selectedFilename = files[selection - 1];
                selectedFilepath = std::filesystem::path(directory) / std::filesystem::path(selectedFilename);
                std::cout << "You selected filepath: " << selectedFilepath << '\n';
            } else {
                std::cout << "\nInvalid selection.\n";
            }
        } else {
            std::cout << "Invalid input. Please enter a valid number.\n";
            std::cin.clear(); // Clear the error state
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Discard invalid input
        }
    } while (!validSelection);

    return selectedFilepath;
}


void readDstr(const std::string& dstrFilepath, std::vector<Particles>& dstr) {
        
    h5::File distribution_file(dstrFilepath, h5::File::ReadOnly);
    //Vectors to save temporarily
    std::vector<real> id_0, latitude_0, alpha_0, aeq_0, ppar_0, pper_0, upar_0, uper_0, Ekin_0, time_0, zeta_0, eta_0, M_adiabatic_0, trapped_0, escaped_0, nan_0, negative_0, high_0;
    //Read dataset from h5file.
    h5::DataSet data_id = distribution_file.getDataSet("id0");
    h5::DataSet data_lat = distribution_file.getDataSet("latitude0");
    h5::DataSet data_aeq = distribution_file.getDataSet("aeq0");
    h5::DataSet data_alpha = distribution_file.getDataSet("alpha0");
    h5::DataSet data_upar = distribution_file.getDataSet("upar0");
    h5::DataSet data_uper = distribution_file.getDataSet("uper0");
    h5::DataSet data_ppar = distribution_file.getDataSet("ppar0");
    h5::DataSet data_pper = distribution_file.getDataSet("pper0");
    h5::DataSet data_eta = distribution_file.getDataSet("eta0");
    h5::DataSet data_zeta = distribution_file.getDataSet("zeta0");
    h5::DataSet data_time = distribution_file.getDataSet("time0");
    h5::DataSet data_M_adiabatic = distribution_file.getDataSet("M_adiabatic0");
    h5::DataSet data_Ekin = distribution_file.getDataSet("Ekin0");
    h5::DataSet data_trapped = distribution_file.getDataSet("trapped0");
    h5::DataSet data_escaped = distribution_file.getDataSet("escaped0");
    h5::DataSet data_nan = distribution_file.getDataSet("nan0");
    h5::DataSet data_negative = distribution_file.getDataSet("negative0");
    h5::DataSet data_high = distribution_file.getDataSet("high0");
    
    // Convert to single vector.
    data_id.read(id_0);
    data_lat.read(latitude_0);
    data_aeq.read(aeq_0);
    data_alpha.read(alpha_0);
    data_upar.read(upar_0);
    data_uper.read(uper_0);
    data_ppar.read(ppar_0);
    data_pper.read(pper_0);
    data_eta.read(eta_0);
    data_zeta.read(zeta_0);
    data_time.read(time_0);
    data_M_adiabatic.read(M_adiabatic_0);
    data_Ekin.read(Ekin_0);
    data_trapped.read(trapped_0);
    data_escaped.read(escaped_0);
    data_escaped.read(nan_0);
    data_escaped.read(negative_0);
    data_escaped.read(high_0);

    // Append to struct from single vector
    for(int p=0; p<Distribution::population; p++)
    {
        dstr[p].id0 = id_0.at(p);
        dstr[p].latitude0 = latitude_0.at(p);
        dstr[p].alpha0 = alpha_0.at(p);  
        dstr[p].aeq0 = aeq_0.at(p);
        dstr[p].ppar0 = ppar_0.at(p);
        dstr[p].pper0 = pper_0.at(p);
        dstr[p].upar0 = upar_0.at(p);
        dstr[p].uper0 = uper_0.at(p);
        dstr[p].Ekin0 = Ekin_0.at(p);
        dstr[p].time0 = time_0.at(p);
        dstr[p].zeta0 = zeta_0.at(p);
        dstr[p].eta0 = eta_0.at(p);
        dstr[p].M_adiabatic0 = M_adiabatic_0.at(p);
        dstr[p].trapped = trapped_0.at(p);
        dstr[p].escaped = escaped_0.at(p);
        dstr[p].nan = nan_0.at(p);
        dstr[p].negative = negative_0.at(p);
        dstr[p].high = high_0.at(p);
    }
}