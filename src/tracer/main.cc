#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <algorithm>
#include <array> 
#include <random> 
#include <iomanip> 
#include <omp.h>
#include <mpi.h>
#include <chrono>
#include <cstdlib> 
#include <ctype.h> 
#include <filesystem>
#include <stdexcept>

// Preprocessor macro instructions are added in files to obey ODR.
#include "no_wpi.h"
#include "bell_wpi.h"
#include "li_wpi.h"
#include "common.h"
#include "struct_Particles.h"   		    	
#include "struct_Ray.h"   		    	
#include "struct_Telescope.h"  				
#include "parameters.h"
#include "constants.h"

//HDF5 I/O
#include <highfive/H5File.hpp>          
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;

int main(int argc, char **argv)
{	
	
	// Initialize MPI
    MPI_Init(&argc, &argv);

    int numProcesses, rank, portionSize;

    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int sizeInBytes[numProcesses];
	int displacements[numProcesses];
	
	// Single particle struct
	Particles single; 

	// Setup arguments
   	SetupArgs setupArgs;

	// Vector of structs for the Particle Distribution
	std::vector<Particles> dstr(Distribution::population, single);	

	// Detector
    Telescope ODPT(Satellite::latitude_deg, Distribution::L_shell);

	// Ray
	Ray ray;
	std::string rayFilepathStr;
	int rayFilepathStrSize;
	
	// Rank 0 reads Particle Distribution and Ray
	if (rank==0) 
	{
	//-------------------------------------------------------- MANAGE COMMAND LINE 	ARGUMENTS  --------------------------------------------------------//

		InputArguments(argc, argv, setupArgs); // If invalid input, terminates program
		
	//------------------------------------------------------------ READ FILES FROM USER  --------------------------------------------------------------//

		std::vector<std::string> dstrFiles, rayFiles;
		std::string selectedFilenameDstr, selectedFilenameRay;
		std::filesystem::path dstrFilepath, rayFilepath;
		const std::string directory = "output/files"; 
		const std::string extension = ".h5"; 
		const std::string dstrPrefix = "dstr"; // Distribution file prefix, defined in distribution.cc
		const std::string rayPrefix = "interpolated_ray"; // Ray file prefix, defined in Read_Ray_Write.cc

		// Read Particle Distribution
		std::cout << "\nPick a Particle Distribution file from the below";
		dstrFilepath = SelectFile(directory, extension, dstrPrefix);
		std::cout << "You selected Particle Distribution filepath: " << dstrFilepath << '\n';
		
		readDstr(dstrFilepath, dstr);

		if (dstr.size()!=Distribution::population) {
			throw std::runtime_error("The parameter for the particle population differs from the population that is defined in the selected hdf5 file! "); 
		}
		else {
			std::cout<<"\nParticle population: "<< Distribution::population <<std::endl;
		}
	
		// If WPI -> Select Ray 
		if(setupArgs.wpi) {
			std::cout << "\nPick a Ray file from the below";
			rayFilepath = SelectFile(directory, extension, rayPrefix);
			std::cout << "You selected Ray filepath: " << rayFilepath << '\n';
			rayFilepathStr = rayFilepath.string(); // Convert to string to broadcast
			rayFilepathStrSize = rayFilepathStr.size();
		}

	//------------------------------------------------------------ DEFINE RANK WORKLOAD  --------------------------------------------------------------//

		// Calculate the size of the portion for each rank
		portionSize = Distribution::population / numProcesses;
		int remainingParticles = Distribution::population % numProcesses;
		// Array to hold the Bytes to be sent to each rank
		sizeInBytes[numProcesses];
		// Array to hold the displacement for each rank
		displacements[numProcesses];
		// Instead of distributing the particle object, we are distributing based on their size in bytes.
		// Distribute in sizeInBytes to distribute univenly with Scatterv
		for (int i = 0; i < numProcesses; ++i) {
			sizeInBytes[i] = portionSize * sizeof(Particles);
		}
		// Distribute remaining particle Bytes evenly 
		for (int i = 0; i < remainingParticles; ++i) {
			sizeInBytes[i] += sizeof(Particles);
		}
		// Calculate the displacement for every rank
		displacements[0] = 0;
		for (int i = 1; i < numProcesses; ++i) {
			displacements[i] = sizeInBytes[i-1] + displacements[i-1];
		}

	}

	//-------------------------------------------------------------- MPI BROADCAST WORKLOAD -----------------------------------------------------------//

	// Master Broadcasts Setup args. (simplified struct broadcast. Works only in single machines/homogeneous clusters)
	MPI_Bcast(&setupArgs, sizeof(SetupArgs), MPI_BYTE, 0, MPI_COMM_WORLD);
	// Master Broadcasts workload data
	MPI_Bcast(&portionSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(sizeInBytes, numProcesses, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displacements, numProcesses, MPI_INT, 0, MPI_COMM_WORLD);
		
	if(setupArgs.wpi) {
		MPI_Bcast(&rayFilepathStrSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		rayFilepathStr.resize(rayFilepathStrSize);
		MPI_Bcast(&rayFilepathStr[0], rayFilepathStrSize, MPI_CHAR, 0, MPI_COMM_WORLD);
		// Construct Ray
		ray.readRay(rayFilepathStr);
	}
	
	//--------------------------------------------------------------- MPI SCATTER WORKLOAD -------------------------------------------------------------//

	dstr.resize(portionSize);

	// Master Scatters distribution unevenly (remaining particles need to be scattered) -> Scatterv
	if (rank==0) {
		MPI_Scatterv(dstr.data(), sizeInBytes, displacements, MPI_BYTE, MPI_IN_PLACE, sizeInBytes[rank], MPI_BYTE, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Scatterv(NULL, sizeInBytes, displacements, MPI_BYTE, dstr.data(), sizeInBytes[rank], MPI_BYTE, 0, MPI_COMM_WORLD);
	}


	for (int i = 0; i < numProcesses; i++) {
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all processes to reach this point
		if (i == rank) {
			std::cout << "Rank " << rank << " particle counts " << sizeInBytes[rank]/sizeof(Particles) << "\n";
			std::cout.flush(); // Ensure that the output is immediately written
		}
	}
	/*
	std::vector<Particles> gathered_dstr;
	if (rank == 0)
	{
		gathered_dstr.resize(Distribution::population);
	}
	MPI_Gatherv(dstr.data(), sizeInBytes[rank], MPI_BYTE, gathered_dstr.data(), sizeInBytes, displacements, MPI_BYTE, 0, MPI_COMM_WORLD);
    // Only Rank 0 prints the gathered data
    if (rank == 0)
    {
        for (const auto& particle : gathered_dstr)
        {
            std::cout << "\ngathered particle " << particle.id0 << " " << particle.latitude0 << std::endl;
        }
    }
	*/

	MPI_Finalize();
	return 0 ;
	
}
	/*
	// Master . Temporary buffer to hold Rank 0's portion of the data
	//std::vector<Particles> dstr_tmp(dstr.begin(), dstr.begin() + portionSize);
	//dstr = dstr_tmp;
	
	// Tests 
	// Ensure all processes have the correct portions of the data
	/*

	}
	// Gather data in Master


//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
	//---NOWPI---//
	if(Nsteps_nowpi>0)
	{

		auto start = std::chrono::high_resolution_clock::now();
		for(int p=displacements[rank]; p<displacements[rank] + sizeInBytes[rank]; p++)     
		{
			// Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
			no_wpi(Nsteps_nowpi, p, dstr[p], ODPT);
		}
		auto stop =  std::chrono::high_resolution_clock::now();

		// Calculate the duration of the execution in seconds (as double)
		std::chrono::duration<real> duration = start - stop;
		std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
		
		// Now initial_particles have the last state after NoWPI
	}
	//---WPI---//
	if(Nsteps_wpi>0)// Run if WPI time is more than 0 seconds
	{

		//---LI---//
		if(argc==3)//Default implementation. Means there's no option '-bell' as command line argument.
		{
			auto start = std::chrono::high_resolution_clock::now();



			for(int p=displacements[rank]; p<displacements[rank] + sizeInBytes[rank]; p++)     
			{
				// Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference
				if(dstr[p].escaped == true) continue; // If this particle is lost, continue with next particle
				if(dstr[p].negative == true) continue; // If this particle came out with aeq negative, continue with next particle
				if(dstr[p].nan == true) continue; // If this particle came out with aeq nan     , continue with next particle
				li_wpi(Nsteps_wpi, p, lat_int, kx_ray, kz_ray, kappa_ray, Bzw, Ezw, Bw_ray, w1, w2, R1, R2, dstr[p], ODPT);  //LI + RAY TRACING
			}
			auto stop =  std::chrono::high_resolution_clock::now();
			
			// Calculate the duration of the execution in seconds (as double)
			std::chrono::duration<real> duration = time2 - time1;
			std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
		}


		//---BELL---//
		else if(argc==4 && argv[3]==bell)
		{
			for(int p=displacements[rank]; p<displacements[rank] + sizeInBytes[rank]; p++)     
			{
				// Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
				if(dstr[p].escaped == true) continue; // If this particle is lost, continue with next particle.
				if(dstr[p].negative == true) continue; // If this particle came out with aeq negative, continue with next particle.
				if(dstr[p].nan == true) continue; // If this particle came out with aeq nan     , continue with next particle.
				bell_wpi(Nsteps_wpi, p, dstr[p], ODPT); // BELL + THE WAVE IS EVERYWHERE
			}
			auto time3 =  std::chrono::high_resolution_clock::now();

			// Calculate the duration of the execution in seconds (as double)
			std::chrono::duration<real> duration = time3 - time2;

			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<duration.count()<<std::endl;

		}
	}


//------------------------------------------------------------ OUTPUT DATA HDF5 --------------------------------------------------------------------------//
	

	// Rank0 saves data to HDF5

	if (rank == 0)
		{ 
		std::cout<<"\nSaving output data to HDF5..."<<std::endl;
		// Assign from struct to vectors.
		std::vector<real> precip_id, precip_latitude, precip_alpha, precip_aeq, precip_time;
		std::vector<real> neg_id, nan_id, high_id;
		std::vector<real> latitude_end, ppar_end, pper_end, alpha_end, aeq_end, eta_end, time_end, Ekin_end;
		std::vector<real> Dsaved_id, Dsaved_max_dEkin, Dsaved_maxEkin_time, Dsaved_max_dPA, Dsaved_maxdPA_time, Dsaved_min_deta_dt, Dsaved_min_deta_dt_time;

		for(int p=0; p<Distribution::population; p++) 
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
		std::string  save_file = "output/files/sim_" + std::to_string(Distribution::population) + "p";  
		if (atoi(argv[1])>0) save_file += std::string("_nowpi") + std::to_string(atoi(argv[1])) + std::string("s"); 
		if (atoi(argv[2])>0) save_file += std::string("_wpi") + std::to_string(atoi(argv[2])) + std::string("s");
		save_file = save_file + ".h5";

		h5::File file(save_file, h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
		
		// Simulation data and Telescope specification - Scalars 
		h5::DataSet telescope_latitude = file.createDataSet("ODPT.latitude_deg", ODPT.latitude_deg);
		h5::DataSet data_population = file.createDataSet("population", Distribution::population);
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

	}


	// Prepare a buffer to hold the received data on rank 0
std::vector<Particles> gatheredData;

// Rank 0 gathers the dstr data from all other ranks
if (rank == 0) {
    // Resize the vector to hold all the data from all ranks
    gatheredData.resize(Distribution::population);
}


	F();


	return 0; 

}
*/
