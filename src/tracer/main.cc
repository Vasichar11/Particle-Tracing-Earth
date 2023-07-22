#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <algorithm>
#include <array> 
#include <random> 
#include <iomanip>  // For std::setprecision()
#include <omp.h>
#include <mpi.h>
#include <chrono>
#include <cstdlib> // atoi()
#include <ctype.h> // isdigit()
#include <filesystem>
#include <stdexcept>

// Preprocessor macro instructions are added in files to obey ODR.
#include "no_wpi.h"
#include "bell_wpi.h"
#include "li_wpi.h"
#include "common.h"
#include "struct_Particles.h"   		    	
#include "struct_Telescope.h"  				
#include "constants.h"

//HDF5 I/O
#include <highfive/H5File.hpp>          
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;

int main(int argc, char **argv)
{	
//------------------------------------------------------------ MANAGE COMMAND LINE 	ARGUMENTS  --------------------------------------------------------------//
	
	// Initialize MPI
    MPI_Init(&argc, &argv);

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// All ranks create a particle distribution
	std::vector<Particles> dstr

	// Only rank 0, the master process, reads the distribution file and constracts the dstr vector
	if (rank==0) 
	{
		
		// Validate number of arguments
		if (argc < 3  || argc > 4) {
			std::cerr << "Invalid number of arguments. Usage: program_name <noWPI_time> <WPI_time> [-bell]\n"
					<< "\n- Argument variables are the simulation times for the NoWPI and WPI interaction in seconds\n"
					<< "- First, provide the time(seconds) of NoWPI simulation and then the time(seconds) of the WPI simulation\n"
					<< "- The default WPI simulation is the one using the Li Formulas (li_wpi.cc)\n"
					<< "- For WPI using the Bell Formulas (bell_wpi.cc), run with the option -bell as the last command line argument\n";
			return 1;
		}


		std::string arg1 = argv[1];
		std::string arg2 = argv[2];

		// Validate arg1 and arg2 as numbers
		try {
			float number1 = std::stof(arg1); // Will throw exception if arguments are not numbers
			float number2 = std::stof(arg2);
			(void)number1;  // Indicate that they are intentionally unused to avoid warnings
			(void)number2;  
		} catch (const std::invalid_argument&) {
			std::cerr << "Invalid argument. Argument 1 and argument 2 must be number of seconds .\n"
					<< "\nArgument variables are the simulation times for the NoWPI and WPI interaction in seconds.\n"
					<< "First, provide the time(seconds) of NoWPI simulation and then the time(seconds) of the WPI simulation.\n"
					<< std::endl;
			return 1;
		}

		std::string bell = "-bell";
		if (argc == 4) {
			std::string arg3 = argv[3];
			if (arg3 != bell) {
				std::cout << "Invalid argument. Argument 3 may only be \"-bell\".\n"
						<< "The default WPI simulation is the one using the Li Formulas (li_wpi.cc)\n"
						<< "For WPI using the Bell Formulas (bell_wpi.cc), run with the option -bell as the last command line argument.\n"
						<< std::endl;
				return 1;
			}
		}
		
		// Assign times from the command line arguments
		real t_nowpi = std::stof(argv[1]); //No WPI time from command line 
		real t_wpi   = std::stof(argv[2]); //WPI time from command line
		real t = t_nowpi + t_wpi;  //Total simulation time
		const int64_t Nsteps_wpi  = t_wpi/Simulation::h; //WPI step count
		const int64_t Nsteps_nowpi= t_nowpi/Simulation::h; //noWPI step count



	//------------------------------------------------------------ READ FILES FROM USER  --------------------------------------------------------------//


		std::cout << "\nPick a particle distribution file from the below\n";
		std::vector<std::string> dstrFiles, rayFiles;
		const std::string directory = "output/files"; // Input file directory
		const std::string extension = ".h5"; // Input file needs to be hdf5 format
		const std::string dstrPrefix = "dstr"; // Input file prefix, defined in distribution.cc
		const std::string rayPrefix = "interpolated_ray"; // Input file prefix, defined in Read_Ray_Write.cc
		std::string selectedFilenameDstr, selectedFilenameRay;
		std::filesystem::path selectedFilepathDstr, selectedFilepathRay;

		// Iterate over the directory and store files that match the distribution hdf5 files
		for (const auto& entry : std::filesystem::directory_iterator(directory)) {
			const std::string filename = entry.path().filename().string();
			const std::string stem = entry.path().stem().string();

			if (entry.path().extension() == extension && stem.substr(0, dstrPrefix.size()) == dstrPrefix) {
				dstrFiles.push_back(filename);
			}
			if (entry.path().extension() == extension && stem.substr(0, rayPrefix.size()) == rayPrefix) {
				rayFiles.push_back(filename);
			}
		}
		// Check if no files were found
		if (dstrFiles.empty()) {
			std::cerr << "Error: No particle distribution files found in \"" << directory << "\".\nPlease create a particle distribution file with distribution.cc and try again." << std::endl;
			return -1;
		}
		int selection;
		bool validSelection1 = false;
		bool validSelection2 = false;

		// Loop until valid particle distribution file selection
		do {
			// Display the files
			std::cout << "\nParticle distribution files in directory \"" << directory << "\":\n";
			for (size_t i = 0; i < dstrFiles.size(); ++i) {
				std::cout << i + 1 << ". " << dstrFiles[i] << '\n';
			}

			// Prompt the user to select a file
			std::cout << "Enter the number of the file you want to select: ";

			if (std::cin >> selection) {
				if (selection >= 1 && selection <= static_cast<int>(dstrFiles.size())) {
					validSelection1 = true;
					selectedFilenameDstr = dstrFiles[selection - 1];
					selectedFilepathDstr = std::filesystem::path(directory) / std::filesystem::path(selectedFilenameDstr);
					std::cout << "You selected filepath: " << selectedFilepathDstr << '\n';
					break;
				} else {
					std::cout << "\nInvalid selection.\n";
				}
			} else {
				std::cout << "Invalid input. Please enter a valid number.\n";
				std::cin.clear(); // Clear the error state
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Discard invalid input
			}
		} while (!validSelection1);

		
		// If WPI then select from the interpolated ray files
		if(t_wpi>0)
		{
		// Loop until valid interpolated ray file selection
			do {
				// Display the files
				std::cout << "\nInterpolated ray files in directory \"" << directory << "\":\n";
				for (size_t i = 0; i < rayFiles.size(); ++i) {
					std::cout << i + 1 << ". " << rayFiles[i] << '\n';
				}

				// Prompt the user to select a file
				std::cout << "Enter the number of the file you want to select: ";

				if (std::cin >> selection) {
					if (selection >= 1 && selection <= static_cast<int>(rayFiles.size())) {
						validSelection2 = true;
						selectedFilenameRay = rayFiles[selection - 1];
						selectedFilepathRay = std::filesystem::path(directory) / std::filesystem::path(selectedFilenameRay);
						std::cout << "You selected filepath: " << selectedFilepathRay << '\n';
						break;
					} else {
						std::cout << "\nInvalid selection.\n";
					}
				} else {
					std::cout << "Invalid input. Please enter a valid number.\n";
					std::cin.clear(); // Clear the error state
					std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Discard invalid input
				}
			} while (!validSelection2);

			std::cout << "You selected filepath: " << selectedFilepathRay << '\n';
		}



	//------------------------------------------------------------READ AND ASSIGN DISTRIBUTION FROM H5 FILE --------------------------------------------------------------//

		h5::File distribution_file(selectedFilepathDstr, h5::File::ReadOnly);
		//Vectors to save temporarily
		std::vector<real> latitude_0, alpha_0, aeq_0, ppar_0, pper_0, upar_0, uper_0, Ekin_0, time_0, zeta_0, eta_0, M_adiabatic_0, trapped_0, escaped_0, nan_0, negative_0, high_0;
		//Read dataset from h5file.
		h5::DataSet data_lat 	     = distribution_file.getDataSet("latitude0");
		h5::DataSet data_aeq	     = distribution_file.getDataSet("aeq0");
		h5::DataSet data_alpha 		 = distribution_file.getDataSet("alpha0");
		h5::DataSet data_upar 		 = distribution_file.getDataSet("upar0");
		h5::DataSet data_uper 		 = distribution_file.getDataSet("uper0");
		h5::DataSet data_ppar 		 = distribution_file.getDataSet("ppar0");
		h5::DataSet data_pper 		 = distribution_file.getDataSet("pper0");
		h5::DataSet data_eta 		 = distribution_file.getDataSet("eta0");
		h5::DataSet data_zeta 		 = distribution_file.getDataSet("zeta0");
		h5::DataSet data_time 		 = distribution_file.getDataSet("time0");
		h5::DataSet data_M_adiabatic = distribution_file.getDataSet("M_adiabatic0");
		h5::DataSet data_Ekin 		 = distribution_file.getDataSet("Ekin0");
		h5::DataSet data_trapped     = distribution_file.getDataSet("trapped0");
		h5::DataSet data_escaped	 = distribution_file.getDataSet("escaped0");
		h5::DataSet data_nan	 	 = distribution_file.getDataSet("nan0");
		h5::DataSet data_negative	 = distribution_file.getDataSet("negative0");
		h5::DataSet data_high    	 = distribution_file.getDataSet("high0");
		//Convert to single vector.
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


		int Population = latitude_0.size(); // Take the population from the h5 file to avoid mistakes (if value in header file changed unitentionally)
		// Single particle struct
		Particles single; 
		// Vector of structs. Has initial states and function to save particle states in vectors
		std::vector<Particles> dstr(Population, single);	

		// Append to struct from single vector
		for(int p=0; p<Population; p++)
		{
			dstr[p].latitude0  	  	  = latitude_0.at(p);
			dstr[p].alpha0  	  	  = alpha_0.at(p);  
			dstr[p].aeq0    	  	  = aeq_0.at(p);
			dstr[p].ppar0   	  	  = ppar_0.at(p);
			dstr[p].pper0   	  	  = pper_0.at(p);
			dstr[p].upar0   	  	  = upar_0.at(p);
			dstr[p].uper0   	  	  = uper_0.at(p);
			dstr[p].Ekin0   	  	  = Ekin_0.at(p);
			dstr[p].time0   	  	  = time_0.at(p);
			dstr[p].zeta0   	  	  = zeta_0.at(p);
			dstr[p].eta0    	  	  = eta_0.at(p);
			dstr[p].M_adiabatic0  	  = M_adiabatic_0.at(p);
			dstr[p].trapped           = trapped_0.at(p);
			dstr[p].escaped           = escaped_0.at(p);
			dstr[p].nan           	  = nan_0.at(p);
			dstr[p].negative          = negative_0.at(p);
			dstr[p].high          	  = high_0.at(p);
		}
		std::cout<<"\nParticle population: "<< Population <<std::endl;


    	// After constructing dstr vector on rank 0

		// Calculate the size of the portion for each rank
		int portionSize = Population / numProcesses;
		
		// Array to hold the number of elements to be sent to each rank
		int counts[numProcesses];
		// Array to hold the displacement for each rank
		int displacements[numProcesses];

		// Populate the counts and displacements arrays based on the portion size
		for (int i = 0; i < numProcesses; ++i) {
			counts[i] = portionSize;
			displacements[i] = i * portionSize;
		}

		// Handle any remaining particles if the population is not divisible evenly
		int remainingParticles = Population % numProcesses;
		for (int i = 0; i < remainingParticles; ++i) {
			counts[i]++;
		}

		// Rank 0 broadcasts the particle population number
	    MPI_Bcast(&Population, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Rank 0 scatters the dstr vector to all other ranks based on the counts and displacements arrays
		MPI_Scatter(dstr.data(), portionSize * sizeof(Particles), MPI_BYTE, dstr.data(), counts[rank] * sizeof(Particles), MPI_BYTE, 0, MPI_COMM_WORLD);
	}
	// All other ranks receive 
	else {
		// Receive the portion of the dstr vector assigned to this rank
   		MPI_Scatter(NULL, 0, MPI_BYTE, dstr.data(), dstr.size() * sizeof(Particles), MPI_BYTE, 0, MPI_COMM_WORLD);
        // Receive the total number of particles from rank 0
        MPI_Bcast(&Population, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // Resize the dstr vector to the appropriate size to hold the portion of data for this process
        dstr.resize(Population);
	}

    // All processes create their local ODPT object
    Telescope ODPT(Satellite::latitude_deg, Distribution::L_shell);

return 0

//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
	int particlesPerRank = Population / numProcesses;
	real wtime = std::chrono::high_resolution_clock::now();

	std::cout<<"Loss cone angle, for Lshell "<<Distribution::L_shell<<" is " << Simulation::alpha_lc*Universal::R2D;

	// ---Selection of THREAD NUM---
	// Reason for being hardcoded
	// Personally I was executing 2 simulations in parallel, both for time "t"
	// That's the reason omp_set_num_threads() is used to hardcode the number of threads
	// e.g Simulation1 "./tracer 10 5" (t=15)
	// e.g Simulation2 "./tracer 15 0" (t=15)
	// Simulation1 for noWPI(8threads) + WPI(20threads) for t = t_noWPI + t_WPI 
	// Simulation2 for noWPI(8threads) for t = t_noWPI
	// This is done to find the effect of the wave particle interaction
	// This selection of threads ensures processors have always work to do
	// After testing the scalability the number of these threads is efficient 
	// Visually
	// Simulation 1
	// start------------|t_noWPI(8threads)|-------end-start-----------------|t_WPI(20threads)|----------------------end
	// Simulation 2
	// start------------|t(8threads)|---------------------------------------------end
	// Overall
	// start------------|16 threads in use| ------end-start--|28 threads in use|--end-start---|20 threads in use|---end

	//---NOWPI---//
		std::cout<<"\n\n"<<t_nowpi<<" sec NoWPI Simulation"<<std::endl;
		std::cout<<"\nExecution time estimation for 8 THREAD run: "<<(Population*0.008/60) * t_nowpi <<" minutes."<<std::endl;
		std::cout<<"Execution time estimation for 20 THREAD run: "<<(Population*0.017/60) * t_nowpi <<" minutes."<<std::endl;

		for(int p=0; p<Population; p++)     // dynamic because some chunks may have less workload.(particles can precipitate and break out of loop)
		{
			// Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
			no_wpi(Nsteps_nowpi, p, dstr[p], ODPT);
		}
		real time1 =  std::chrono::high_resolution_clock::now();
		std::cout<<"\nExecution time: "<<time1<<std::endl;

	// Now initial_particles have the last state after NoWPI

	//---WPI---//
	if(t_wpi>0)// Run if WPI time is more than 0 seconds
	{

		//---LI---//
		if(argc==3)//Default implementation. Means there's no option '-bell' as command line argument.
		{
			std::cout<<"\n\n"<<t_wpi<<" sec NoWPI Simulation (Li+ray tracing)."<<std::endl;
			// Read hdf5 files from disk to pass them to the WPI function instead of reading them in every thread - every time - accessing disk frequently
			// read_hdf5() is a function to read HDF5 dataset as vectors. 
			
			//std::cout << "\nPick Ray file from the below\n";
			//for (const auto & entry : std::filesystem::directory_iterator("output/files")) {
			//	std::cout << entry.path() << std::endl; }
			std::string file_ray;
			std::cout<<"\n";
			//std::getline(std::cin, file_ray);
			file_ray = "output/files/interpolated_ray_pwr1.000000.h5";
			std::cout << "\n\nRay file: "<<file_ray;
			

			std::vector <real> lat_int = read_hdf5("lat_int", file_ray);
			const std::vector <real> kx_ray = read_hdf5("kx_ray", file_ray);    
			const std::vector <real> kz_ray = read_hdf5("kz_ray", file_ray);   
			const std::vector <real> kappa_ray = read_hdf5("kappa_ray", file_ray);       
			const std::vector <real> Bzw = read_hdf5("Bzw", file_ray);
			const std::vector <real> Ezw = read_hdf5("Ezw", file_ray);
			const std::vector <real> Bw_ray = read_hdf5("Bw_ray", file_ray);    
			const std::vector <real> w1 = read_hdf5("w1", file_ray);
			const std::vector <real> w2 = read_hdf5("w2", file_ray);
			const std::vector <real> R1 = read_hdf5("R1", file_ray);
			const std::vector <real> R2 = read_hdf5("R2", file_ray);

			for(int p=0; p<Population; p++)     
			{
				// Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference
				if(dstr[p].escaped  == true) continue; // If this particle is lost, continue with next particle
				if(dstr[p].negative == true) continue; // If this particle came out with aeq negative, continue with next particle
				if(dstr[p].nan 	    == true) continue; // If this particle came out with aeq nan     , continue with next particle
				li_wpi(Nsteps_wpi, p, lat_int, kx_ray, kz_ray, kappa_ray, Bzw, Ezw, Bw_ray, w1, w2, R1, R2, dstr[p], ODPT);  //LI   + RAY TRACING
			}
			real time2 = std::chrono::high_resolution_clock::now();
			std::cout<<"\nExecution time: "<<time2<<std::endl;
		}

		//---BELL---//
		else if(argc==4 && argv[3]==bell)
		{
			std::cout<<"\n\n"<<t_wpi<<" sec NoWPI Simulation(Bell).\nWave magnitude(T): "<<Wave::By_wave<<std::endl;
			std::cout<<"Execution time estimation for 8 THREAD run: "<<(Population*0.036/60) * t_wpi <<" minutes."<<std::endl;
			for(int p=0; p<Population; p++)     
			{
				// Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
				if(dstr[p].escaped  == true) continue; // If this particle is lost, continue with next particle.
				if(dstr[p].negative == true) continue; // If this particle came out with aeq negative, continue with next particle.
				if(dstr[p].nan 	    == true) continue; // If this particle came out with aeq nan     , continue with next particle.
				bell_wpi(Nsteps_wpi, p, dstr[p], ODPT); 		// BELL + THE WAVE IS EVERYWHERE
			}
			real time2 =  std::chrono::high_resolution_clock::now();
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}
	}


//------------------------------------------------------------ OUTPUT DATA HDF5 --------------------------------------------------------------------------//
	std::cout<<"\nSaving output data to HDF5..."<<std::endl;
	// Assign from struct to vectors.
	std::vector<real> precip_id, precip_latitude, precip_alpha, precip_aeq, precip_time;
	std::vector<real> neg_id, nan_id, high_id;
	std::vector<real> latitude_end, ppar_end, pper_end, alpha_end, aeq_end, eta_end, time_end, Ekin_end;
	std::vector<real> Dsaved_id, Dsaved_max_dEkin, Dsaved_maxEkin_time, Dsaved_max_dPA, Dsaved_maxdPA_time, Dsaved_min_deta_dt, Dsaved_min_deta_dt_time;

	for(int p=0; p<Population; p++) 
	{
		// Last particle states(that can become first states for next simulation).
		latitude_end.push_back(dstr[p].latitude_end);
    	ppar_end.push_back(dstr[p].ppar_end); 
    	pper_end.push_back(dstr[p].pper_end); 
    	alpha_end.push_back(dstr[p].alpha_end); 
    	aeq_end.push_back(dstr[p].aeq_end); 
    	eta_end.push_back(dstr[p].eta_end); 
    	Ekin_end.push_back(dstr[p].Ekin_end); 
    	time_end.push_back(dstr[p].time_end);
		

		// Precipitating Particles
		if(dstr[p].escaped) 
		{
			precip_id.push_back(dstr[p].id_lost);
			precip_latitude.push_back(dstr[p].latitude_lost); 
			precip_alpha.push_back(dstr[p].alpha_lost);
			precip_aeq.push_back(dstr[p].aeq_lost);
			precip_time.push_back(dstr[p].time_lost);
		}

		// Negative P.A Particles
		if(dstr[p].negative) 
		{
			neg_id.push_back(dstr[p].id_neg);
		}
		// High P.A Particles
		if(dstr[p].high) 
		{
			high_id.push_back(dstr[p].id_neg);
		}
		// High P.A Particles
		if(dstr[p].nan) 
		{
			nan_id.push_back(dstr[p].id_nan);
		}

		// Saved states 
		Dsaved_id.push_back(dstr[p].saved_id); 
		Dsaved_max_dEkin.push_back(dstr[p].saved_max_dEkin); 
		Dsaved_maxEkin_time.push_back(dstr[p].saved_maxEkin_time); 
		Dsaved_max_dPA.push_back(dstr[p].saved_max_dPA); 
		Dsaved_maxdPA_time.push_back(dstr[p].saved_maxdPA_time); 
		Dsaved_min_deta_dt.push_back(dstr[p].saved_min_detadt); 
		Dsaved_min_deta_dt_time.push_back(dstr[p].saved_mindetadt_time); 
	}


	// File name based on the argument variables
	std::string  save_file = "output/files/sim_" + std::to_string(Population) + "p";  
	if (atoi(argv[1])>0) save_file += std::string("_nowpi") + std::to_string(atoi(argv[1])) + std::string("s"); 
	if (atoi(argv[2])>0) save_file += std::string("_wpi") + std::to_string(atoi(argv[2])) + std::string("s");
	save_file = save_file + ".h5";

	h5::File file(save_file, h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	// Simulation data and Telescope specification - Scalars 
	h5::DataSet telescope_latitude = file.createDataSet("ODPT.latitude_deg", ODPT.latitude_deg);
	h5::DataSet data_population = file.createDataSet("population", Population);
	h5::DataSet simulation_time = file.createDataSet("t", t);
	
	// Detected particles from the satellite
	h5::DataSet detected_latitude = file.createDataSet("ODPT.latitude", ODPT.latitude);
	h5::DataSet detected_time = file.createDataSet("ODPT.time", ODPT.time);
	h5::DataSet detected_id = file.createDataSet("ODPT.id", ODPT.id);
	h5::DataSet detected_alpha = file.createDataSet("ODPT.alpha", ODPT.alpha);
	h5::DataSet detected_aeq = file.createDataSet("ODPT.aeq", ODPT.aeq); 

	// Precipitating Particles
	h5::DataSet precipitated_id = file.createDataSet("precip_id", precip_id);
	h5::DataSet precipitated_latitude = file.createDataSet("precip_latitude", precip_latitude);
	h5::DataSet precipitated_aeq = file.createDataSet("precip_aeq", precip_aeq);
	h5::DataSet precipitated_alpha = file.createDataSet("precip_alpha", precip_alpha);
	h5::DataSet precipitated_time = file.createDataSet("precip_time", precip_time);

	// Negative P.A Particles
	h5::DataSet data_neg_id = file.createDataSet("neg_id", neg_id);
	// High P.A Particles
	h5::DataSet data_high_id = file.createDataSet("high_id", high_id);
	// High P.A Particles
	h5::DataSet data_nan_id = file.createDataSet("nan_id", nan_id);

	// Particles states at noWPI end.
	h5::DataSet ending_latutide = file.createDataSet("latitude_end", latitude_end);
	h5::DataSet ending_ppar = file.createDataSet("ppar_end", ppar_end);
	h5::DataSet ending_pper = file.createDataSet("pper_end", pper_end);
	h5::DataSet ending_alpha = file.createDataSet("alpha_end", alpha_end);
	h5::DataSet ending_aeq = file.createDataSet("aeq_end", aeq_end);
	h5::DataSet ending_eta = file.createDataSet("eta_end", eta_end);
	h5::DataSet ending_Ekin = file.createDataSet("Ekin_end", Ekin_end);
	h5::DataSet ending_time = file.createDataSet("time_end", time_end);

	// Saved Particles
	h5::DataSet saved_id = file.createDataSet("saved_id", Dsaved_id);
	h5::DataSet saved_max_dEkin = file.createDataSet("saved_max_dEkin", Dsaved_max_dEkin);
	h5::DataSet saved_maxEkin_time = file.createDataSet("saved_maxEkin_time", Dsaved_maxEkin_time);
	h5::DataSet saved_max_dPA = file.createDataSet("saved_max_dPA", Dsaved_max_dPA);
	h5::DataSet saved_maxdPA_time = file.createDataSet("saved_maxdPA_time", Dsaved_maxdPA_time);
	h5::DataSet saved_mindeta_dt = file.createDataSet("saved_mindeta_dt", Dsaved_min_deta_dt);
	h5::DataSet saved_mindeta_dt_time = file.createDataSet("saved_mindeta_dt_time", Dsaved_min_deta_dt_time);
	std::cout<<"\nData saved!"<<std::endl;

//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//
	std::cout<<"\nSimulation end"<<std::endl;

return 0; 

}


