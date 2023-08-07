#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <algorithm>
#include <array> 
#include <random> 
#include <iomanip> 
#include <omp.h>
#include <chrono>
#include <cstdlib> 
#include <ctype.h> 
#include <filesystem>
#include <stdexcept>

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
	// Single particle struct
	Particles particle; 
	
	// Vector of structs for the Particle Distribution
	std::vector<Particles> dstr(Distribution::population, particle);	

	// Setup arguments
   	SetupArgs setupArgs;

	// Detector
	Telescope ODPT(Satellite::latitude_deg, Distribution::L_shell);
	
	// Ray
	Ray ray;

	//-------------------------------------------------------- MANAGE COMMAND LINE 	ARGUMENTS  --------------------------------------------------------//

		InputArguments(argc, argv, setupArgs); // If invalid input, terminates program
		
	//------------------------------------------------------------ READ FILES FROM USER  --------------------------------------------------------------//

		std::vector<std::string> dstrFiles, rayFiles;
		std::string selectedFilenameDstr, selectedFilenameRay, rayFilepathStr;
		std::filesystem::path dstrFilepath, rayFilepath;
		std::filesystem::path directory = std::filesystem::path("output") / "files";
		const std::string extension = ".h5"; 
		const std::string dstrPrefix = "dstr"; // Distribution file prefix, defined in distribution.cc
		const std::string rayPrefix = "interpolated_ray"; // Ray file prefix, defined in Read_Ray_Write.cc

		// Read Particle Distribution
		std::cout << "\nPick a Particle Distribution file:\n";
		dstrFilepath = SelectFile(directory, extension, dstrPrefix);
		
		readDstr(dstrFilepath, dstr);

		if (dstr.size()!=Distribution::population) {
			throw std::runtime_error("The parameter for the particle population differs from the population that is defined in the selected hdf5 file! "); 
		}
		else {
			std::cout<<"\nParticle population: "<< Distribution::population <<std::endl;
		}
	
		// If WPI -> Select Ray 
		if(setupArgs.wpi) {
			std::cout << "\nPick a Ray file:\n";
			rayFilepath = SelectFile(directory, extension, rayPrefix);
			rayFilepathStr = rayFilepath.string(); // Convert to string to broadcast
			ray.readRay(rayFilepathStr);
		}
//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
    std::cout << "\nSimulation starts...\n";
	int realthreads;   
	real wtime = omp_get_wtime();
	
	//---NOWPI---//
	if(setupArgs.nowpi) {
		
		omp_set_num_threads(Simulation::num_threads_nowpi); // Performance peaks with 8 threads. 
		std::cout<<"\n\n"<<setupArgs.t_nowpi<<" sec NoWPI Simulation"<<std::endl;
		std::cout<<"\nExecution time estimation for 8 THREAD run: "<<(Distribution::population*0.008/60) * setupArgs.t_nowpi <<" minutes."<<std::endl;
		std::cout<<"Execution time estimation for 20 THREAD run: "<<(Distribution::population*0.017/60) * setupArgs.t_nowpi <<" minutes."<<std::endl;
		std::cout<<"\nForked...";

		#pragma omp parallel
		{
			int id = omp_get_thread_num();
			if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
			#pragma omp for schedule(dynamic)
			for (auto &particle : dstr) {
				int id = std::distance(dstr.begin(), std::find(dstr.begin(), dstr.end(), particle));
				no_wpi(setupArgs.Nsteps_nowpi, id, particle, ODPT);
			}
		}
		
		std::cout<<"\n\n"<<"Joined"<<std::endl;
		real time1 = omp_get_wtime()-wtime;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time1<<std::endl;
	}

	//---WPI---//
	if(setupArgs.wpi) {
		
		//---LI---//
		if(setupArgs.use_li_equations) {

			omp_set_num_threads(Simulation::num_threads_wpi_li); // Better performance for more threads
			std::cout<<"\n\n"<<setupArgs.t_wpi<<" sec NoWPI Simulation (Li+ray tracing)."<<std::endl;
			std::cout<<"Execution time estimation for 20 THREAD run: "<<(Distribution::population*0.1/60) * setupArgs.t_wpi <<" minutes."<<std::endl;
			std::cout<<"\nForked...";

			#pragma omp parallel
			{
 				int id = omp_get_thread_num();
				if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
				#pragma omp for schedule(dynamic)
				for (auto &particle : dstr) {
					if(particle.escaped || particle.negative || particle.nan) continue; // If particle was lost, negative or nan, continue with next particle 
					int id = std::distance(dstr.begin(), std::find(dstr.begin(), dstr.end(), particle));
					li_wpi(setupArgs.Nsteps_wpi, id, particle, ODPT, ray);  
				}
			}

			std::cout<<"\n\n"<<"Joined"<<std::endl;
			real time2 = omp_get_wtime()-wtime ;
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}

		//---BELL---//
		else if(setupArgs.use_bell_equations) {

			omp_set_num_threads(Simulation::num_threads_wpi_bell); // Better performance for more threads
			std::cout<<"\n\n"<<setupArgs.t_wpi<<" sec NoWPI Simulation(Bell).\nWave magnitude(T): "<<Wave::By_wave<<std::endl;
			std::cout<<"Execution time estimation for 8 THREAD run: "<<(Distribution::population*0.036/60) * setupArgs.t_wpi <<" minutes."<<std::endl;
			std::cout<<"\nForked...";

			#pragma omp parallel
			{
				int id = omp_get_thread_num();
				if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
				#pragma omp for schedule(dynamic)
				for (auto &particle : dstr) {
					if(particle.escaped || particle.negative || particle.nan) continue; // If particle was lost, negative or nan, continue with next particle 
					int id = std::distance(dstr.begin(), std::find(dstr.begin(), dstr.end(), particle));
					bell_wpi(setupArgs.Nsteps_wpi, id, particle, ODPT);  
				}
			}

			std::cout<<"\n\n"<<"Joined"<<std::endl;
			real time2 = omp_get_wtime()-wtime ;
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}
	}

//-------------------------------------------------------------------- SAVE DATA --------------------------------------------------------------------------//
	
	
	std::cout<<"\nStrings: "<<std::to_string(setupArgs.t_nowpi)<< " " << std::to_string(setupArgs.t_wpi)<<std::endl;
	std::cout<<"\nSaving output data to HDF5..."<<std::endl;
	// Assign from struct to vectors.
	std::vector<real> precip_id, precip_latitude, precip_alpha, precip_aeq, precip_time;
	std::vector<real> neg_id, nan_id, high_id;
	std::vector<real> latitude_end, ppar_end, pper_end, alpha_end, aeq_end, eta_end, time_end, Ekin_end;
	std::vector<real> Dsaved_id, Dsaved_max_dEkin, Dsaved_maxEkin_time, Dsaved_max_dPA, Dsaved_maxdPA_time, Dsaved_min_deta_dt, Dsaved_min_deta_dt_time;

	for(int p=0; p<Distribution::population; p++) {
		
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
		if(dstr[p].escaped) {
			precip_id.push_back(dstr[p].id_lost);
			precip_latitude.push_back(dstr[p].latitude_lost); 
			precip_alpha.push_back(dstr[p].alpha_lost);
			precip_aeq.push_back(dstr[p].aeq_lost);
			precip_time.push_back(dstr[p].time_lost);
		}
		
		// Negative P.A Particles
		if(dstr[p].negative) neg_id.push_back(dstr[p].id_neg);
		// High P.A Particles
		if(dstr[p].high) high_id.push_back(dstr[p].id_neg);
		// High P.A Particles
		if(dstr[p].nan) nan_id.push_back(dstr[p].id_nan);

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
	std::filesystem::path outputFile = directory / ("sim_" + std::to_string(Distribution::population) + "p");
	
	std::ostringstream t_nowpi, t_wpi;
    t_nowpi << setupArgs.t_nowpi; 
    t_wpi << setupArgs.t_wpi; 

	if (setupArgs.nowpi) outputFile += std::string("_nowpi") + t_nowpi.str()+ std::string("s");
	if (setupArgs.wpi) outputFile += std::string("_wpi") + t_wpi.str() + std::string("s");

	outputFile += std::string(".h5");

	// To get the string representation of the path
	std::string outputFile_str = outputFile.string();

	h5::File file(outputFile, h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	// Simulation data and Telescope specification - Scalars 
	h5::DataSet telescope_latitude = file.createDataSet("ODPT.latitude_deg", ODPT.latitude_deg);
	h5::DataSet data_population = file.createDataSet("population", Distribution::population);
	h5::DataSet simulation_time = file.createDataSet("t", setupArgs.t);
	
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
	std::cout<<"\nSimulation end"<<std::endl;

	return 0; 
	
}