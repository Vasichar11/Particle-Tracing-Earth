#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <algorithm>
#include <array> 
#include <random> 
#include <iomanip>  //For std::setprecision()
#include <omp.h>
#include <cstdlib> //atoi()
#include <ctype.h> //isdigit()

//Same directory headers							    
//Preprocessor macro instructions are added in files to obey ODR.
#include "headers/no_wpi.h"
#include "headers/bell_wpi.h"
#include "headers/li_wpi.h"
#include <filesystem>

#include "headers/common.h"
#include "headers/struct_Particles.h"   		    	
#include "headers/struct_Telescope.h"  				
#include "headers/constants.h"

//HDF5 I/O
#include <highfive/H5File.hpp>          
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;


int main(int argc, char **argv)
{
	//Simulation time and steps, read from command line.
	real t_nowpi = atoi(argv[1]);    //No WPI time from command line 
	real t_wpi   = atoi(argv[2]);	 //WPI time from command line
	real t = t_nowpi + t_wpi;  //Total simulation time
	int64_t Nsteps_wpi  = t_wpi/Constants::h; 	 //WPI step count
	int64_t Nsteps_nowpi= t_nowpi/Constants::h; //noWPI step count
	
	std::string bell = "-bell";
	/*//---ARGV ERROR---//
	if( !(isdigit(argv[1])) || !(isdigit(argv[2])) )
	{ 
		std::cout<<"\nArgument variables are the simulation times for the NoWPI and WPI interaction. First provide the time(s) of NoWPI simulation and then the time(s) of the WPI simulation. For WPI using the Bell Formulas run with the option -bell as the last command line argument.\n"<<std::endl;
		return EXIT_FAILURE;	
	}*/
	if(argc==4)
	{
		if(argv[3]!=bell){
			std::cout<<"\nThird argument can be used only to run the simulation using the Bell Formulas, providing the option '-bell' as the last argument\n"<<std::endl;
		}
	}



	std::cout << "\nPick a particle distribution file from the below\n";
	for (const auto & entry : std::filesystem::directory_iterator("h5files")) {
        std::cout << entry.path() << std::endl; }
	std::string str;
	std::cout<<"\n";
	std::getline(std::cin, str);
	std::cout << "\n\nParticle distribution file: "<<str;

//------------------------------------------------------------READ AND ASSIGN DISTRIBUTION FROM H5 FILE --------------------------------------------------------------//
	h5::File distribution_file(str, h5::File::ReadOnly);
	//Vectors to save temporarily
	std::vector<real> lamda_0, alpha_0, aeq_0, ppar_0, pper_0, upar_0, uper_0, Ekin_0, time_0, zeta_0, eta_0, M_adiabatic_0, trapped_0, escaped_0, nan_0, negative_0, high_0;
	//Read dataset from h5file.
	h5::DataSet data_lat 	     = distribution_file.getDataSet("lamda0");
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
	data_lat.read(lamda_0);
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


	int Population = lamda_0.size(); //Take the population from the h5 file to avoid mistakes.
	//Single particle struct.
	Particles single; 
	//Vector of structs. Has initial states and function to save particle states in vectors.
	std::vector<Particles> dstr(Population, single);	

	//Append to struct from single vector.
	for(int p=0; p<Population; p++)
	{
		dstr[p].lamda0  	  	  = lamda_0.at(p);
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


	//Object for Particle Telescope.		
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell);	


//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
	int realthreads;   
	real wtime = omp_get_wtime();


	//---NOWPI---//
	omp_set_num_threads(8); //Performance peaks with 8 threads.
	if(t_nowpi>0) //Run if noWPI time is more than 0 seconds.
	{
		std::cout<<"\n\n"<<t_nowpi<<" sec NoWPI Simulation"<<std::endl;
		std::cout<<"\nExecution time estimation for 8 THREAD run: "<<(Population*0.008/60) * t_nowpi <<" minutes."<<std::endl;
		std::cout<<"Execution time estimation for 20 THREAD run: "<<(Population*0.017/60) * t_nowpi <<" minutes."<<std::endl;
		std::cout<<"\nForked...";
		//---PARALLELISM Work sharing---//
		
		#pragma omp parallel
		{
			int id = omp_get_thread_num();
			if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
			#pragma omp for schedule(dynamic)
				for(int p=0; p<Population; p++)     //dynamic because some chunks may have less workload.(particles can precipitate and break out of loop)
				{
					//Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
					no_wpi(Nsteps_nowpi, p, dstr[p], ODPT);
				}
		}	
		std::cout<<"\n\n"<<"Joined"<<std::endl;
		real time1 = omp_get_wtime()-wtime;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time1<<std::endl;
	}

	//Now initial_particles have the last state after NoWPI

	//---WPI---//
	omp_set_num_threads(20); //Better performance for more threads
	if(t_wpi>0)//Run if WPI time is more than 0 seconds.
	{
		//---LI---//
		if(argc==3)//Default implementation. Means there's no option '-bell' as command line argument.
		{
			std::cout<<"\n\n"<<t_wpi<<" sec ray tracing WPI Simulation using Li formulas."<<std::endl;
			std::cout<<"Execution time estimation for 20 THREAD run: "<<(Population*0.1/60) * t_wpi <<" minutes."<<std::endl;
			std::cout<<"\nForked...";
			//--READ HDF5 files from disk to pass them to the WPI function instead of reading them in every thread - every time - accessing disk frequently --//
			//read_hdf5() is a function to read HDF5 dataset as vectors. 
			std::vector <real> lat_int       =   read_hdf5("lat_int",       "h5files/interpolated_ray.h5");
			const std::vector <real> kx_ray        =   read_hdf5("kx_ray",        "h5files/interpolated_ray.h5");    
			const std::vector <real> kz_ray        =   read_hdf5("kz_ray",        "h5files/interpolated_ray.h5");   
			const std::vector <real> kappa_ray     =   read_hdf5("kappa_ray",     "h5files/interpolated_ray.h5");       
			const std::vector <real> Bzw           =   read_hdf5("Bzw",           "h5files/interpolated_ray.h5");
			const std::vector <real> Ezw           =   read_hdf5("Ezw",           "h5files/interpolated_ray.h5");
			const std::vector <real> Bw_ray        =   read_hdf5("Bw_ray",        "h5files/interpolated_ray.h5");    
			const std::vector <real> w1            =   read_hdf5("w1",            "h5files/interpolated_ray.h5");
			const std::vector <real> w2            =   read_hdf5("w2",            "h5files/interpolated_ray.h5");
			const std::vector <real> R1            =   read_hdf5("R1",            "h5files/interpolated_ray.h5");
			const std::vector <real> R2            =   read_hdf5("R2",            "h5files/interpolated_ray.h5");
			//---PARALLELISM Work sharing---//
			#pragma omp parallel
			{
 				int id = omp_get_thread_num();
				if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
				#pragma omp for schedule(dynamic)
					for(int p=0; p<Population; p++)     
					{
						//Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
						if(dstr[p].escaped  == true) continue; //If this particle is lost, continue with next particle.
						if(dstr[p].negative == true) continue; //If this particle came out with aeq negative, continue with next particle.
						if(dstr[p].nan 	    == true) continue; //If this particle came out with aeq nan     , continue with next particle.
						li_wpi(Nsteps_wpi, p, lat_int, kx_ray, kz_ray, kappa_ray, Bzw, Ezw, Bw_ray, w1, w2, R1, R2, dstr[p], ODPT);  //LI   + RAY TRACING
					}
			}	
			std::cout<<"\n\n"<<"Joined"<<std::endl;
			real time2 = omp_get_wtime()-wtime ;
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}

		//---BELL---//
		else if(argc==4 && argv[3]==bell)
		{
			std::cout<<"\n\n"<<t_wpi<<" sec WPI Simulation using Bell formulas. Wave magnitude(T): "<<Constants::By_wave<<std::endl;
			std::cout<<"Execution time estimation for 8 THREAD run: "<<(Population*0.036/60) * t_wpi <<" minutes."<<std::endl;
			std::cout<<"\nForked...";
			//---PARALLELISM Work sharing---//
			#pragma omp parallel
			{
				int id = omp_get_thread_num();
				if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
				#pragma omp for schedule(dynamic)
					for(int p=0; p<Population; p++)     
					{
						//Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
						if(dstr[p].escaped  == true) continue; //If this particle is lost, continue with next particle.
						if(dstr[p].negative == true) continue; //If this particle came out with aeq negative, continue with next particle.
						if(dstr[p].nan 	    == true) continue; //If this particle came out with aeq nan     , continue with next particle.
						bell_wpi(Nsteps_wpi, p, dstr[p], ODPT); 		//BELL + THE WAVE IS EVERYWHERE
					}
			}	
			std::cout<<"\n\n"<<"Joined"<<std::endl;
			real time2 = omp_get_wtime()-wtime ;
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}
	}
//------------------------------------------------------------------ SIMULATION: END ---------------------------------------------------------------------//

//------------------------------------------------------------ OUTPUT DATA HDF5 --------------------------------------------------------------------------//
 
	//Assign from struct to vectors.
	std::vector<real> precip_id, precip_lamda, precip_alpha, precip_aeq, precip_time, neg_id, nan_id, high_id, lamda00, ppar00, pper00, alpha00, aeq00, eta00, time00, Ekin00;
	std::vector<real> saved_deta_dt,saved_id,saved_lamda, saved_Ekin, saved_alpha, saved_eta;

	for(int p=0; p<Population; p++) 
	{
		//Last particle states(that can become first states for next simulation).
		lamda00.push_back(dstr[p].lamda00);
    	ppar00.push_back(dstr[p].ppar00); 
    	pper00.push_back(dstr[p].pper00); 
    	alpha00.push_back(dstr[p].alpha00); 
    	aeq00.push_back(dstr[p].aeq00); 
    	eta00.push_back(dstr[p].eta00); 
    	Ekin00.push_back(dstr[p].Ekin00); 
    	time00.push_back(dstr[p].time00);
		

		//Precipitating Particles
		if(dstr[p].escaped) 
		{
			precip_id.push_back(dstr[p].id_lost);
			precip_lamda.push_back(dstr[p].lamda_lost); 
			precip_alpha.push_back(dstr[p].alpha_lost);
			precip_aeq.push_back(dstr[p].aeq_lost);
			precip_time.push_back(dstr[p].time_lost);
		}

		//Negative P.A Particles
		if(dstr[p].negative) 
		{
			neg_id.push_back(dstr[p].id_neg);
		}
		//High P.A Particles
		if(dstr[p].high) 
		{
			high_id.push_back(dstr[p].id_neg);
		}
		//High P.A Particles
		if(dstr[p].nan) 
		{
			nan_id.push_back(dstr[p].id_nan);
		}

		//Minimum deta_dt state 
 		saved_deta_dt.push_back(dstr[p].saved_deta_dt);
 		saved_id.push_back(dstr[p].saved_id);
 		saved_lamda.push_back(dstr[p].saved_lamda);
 		saved_Ekin.push_back(dstr[p].saved_Ekin);
 		saved_eta.push_back(dstr[p].saved_eta);
 		saved_alpha.push_back(dstr[p].saved_alpha);

	}


	//File name based on the argument variables
	std::string  save_file = "h5files/" + std::to_string(Constants::population) + "p_";   
	if 	  	(argc==2) save_file = save_file  + std::string(argv[1]) ;
	else if (argc==3) save_file = save_file  + std::string("both"); 
	save_file = save_file + ".h5";

	h5::File file(save_file, h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	//Simulation data and Telescope specification - Scalars 
	h5::DataSet telescope_lamda = file.createDataSet("ODPT.latitude", ODPT.latitude);
	h5::DataSet data_population = file.createDataSet("population", 	Population);
	h5::DataSet simulation_time = file.createDataSet("t", 			t);
	
	//Detected particles from the satellite
	h5::DataSet detected_lamda  = file.createDataSet("ODPT.lamda", ODPT.lamda);
	h5::DataSet detected_time   = file.createDataSet("ODPT.time", ODPT.time);
	h5::DataSet detected_id     = file.createDataSet("ODPT.id", ODPT.id);
	h5::DataSet detected_alpha  = file.createDataSet("ODPT.alpha", ODPT.alpha);
	h5::DataSet detected_aeq    = file.createDataSet("ODPT.aeq", ODPT.aeq); 

	//Precipitating Particles
	h5::DataSet precipitated_id    = file.createDataSet("precip_id", 	 precip_id);
	h5::DataSet precipitated_lamda = file.createDataSet("precip_lamda",  precip_lamda);
	h5::DataSet precipitated_aeq   = file.createDataSet("precip_aeq", 	 precip_aeq);
	h5::DataSet precipitated_alpha = file.createDataSet("precip_alpha",  precip_alpha);
	h5::DataSet precipitated_time  = file.createDataSet("precip_time",	 precip_time);

	//Negative P.A Particles
	h5::DataSet data_neg_id		= file.createDataSet("neg_id", 	 neg_id);
	//High P.A Particles
	h5::DataSet data_high_id	= file.createDataSet("high_id",  high_id);
	//High P.A Particles
	h5::DataSet data_nan_id		= file.createDataSet("nan_id", 	 nan_id);


	//Particles states at noWPI end.
	h5::DataSet ending_lamda   = file.createDataSet("lamda00", lamda00);
	h5::DataSet ending_ppar    = file.createDataSet("ppar00",  ppar00);
	h5::DataSet ending_pper    = file.createDataSet("pper00",  pper00);
	h5::DataSet ending_alpha   = file.createDataSet("alpha00", alpha00);
	h5::DataSet ending_aeq     = file.createDataSet("aeq00",   aeq00);
	h5::DataSet ending_eta     = file.createDataSet("eta00",   eta00);
	h5::DataSet ending_Ekin    = file.createDataSet("Ekin00",  Ekin00);
	h5::DataSet ending_time    = file.createDataSet("time00",  time00);


	//Saved Particles
	h5::DataSet saved_id_data        = file.createDataSet("saved_id", 	  saved_id);
	h5::DataSet saved_lamda_data     = file.createDataSet("saved_lamda",  saved_lamda);
	h5::DataSet saved_deta_dt_data   = file.createDataSet("saved_deta_dt",saved_deta_dt);
	h5::DataSet saved_Ekin_data      = file.createDataSet("saved_Ekin",   saved_Ekin);
	h5::DataSet saved_eta_data       = file.createDataSet("saved_eta",    saved_eta);
	h5::DataSet saved_alpha_data     = file.createDataSet("saved_alpha",  saved_alpha);
//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//
return 0; 

}


