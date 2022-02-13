#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <algorithm>
#include <array> 
#include <random> 
#include <iomanip>  //For std::setprecision()
#include <omp.h>
 
//Same directory headers							    
//Preprocessor macro instructions are added in files to obey ODR.
#include "headers/no_wpi.h"
#include "headers/bell_wpi.h"
#include "headers/li_wpi.h"


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
	std::string string_no_wpi   = "no_wpi";
	std::string string_li_wpi   = "li_wpi";
	std::string string_bell_wpi = "bell_wpi";
	
	//---ARGV ERROR---//
	if(  argv[1]!=string_no_wpi    ||    ( argc==3 && argv[2]!=string_bell_wpi && argv[2]!=string_li_wpi ) )
	{ 
		std::cout<<"\n\nArgument variables don't match any of the program's possible implementations.\nSet first argv=no_wpi and third argv from the list:(bell_wpi, li_wpi).\n"<<std::endl;
		return EXIT_FAILURE;	
	}

	
	//Object for Particle Telescope.		
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell);		
	
	//Single particle struct.
	Particles single; 

	//Vector of structs to store desired particle states.
	std::vector<Particles> lost_particles(Constants::population, single);	

	//Vector of structs for initial particle states.
	std::vector<Particles> initial_particles(Constants::population, single);	

//------------------------------------------------------------READ DISTRIBUTION FROM H5 FILE --------------------------------------------------------------//
	h5::File distribution_file("h5files/distribution_test.h5", h5::File::ReadOnly);
	//Vectors to save temporarily
	std::vector<real> lamda_0, alpha_0, aeq_0, ppar_0, pper_0, upar_0, uper_0, Ekin_0, time_0, zeta_0, eta_0, deta_dt_0, M_adiabatic_0;
	//Read dataset from h5file.
	h5::DataSet data_lat 	     = distribution_file.getDataSet("lat");
	h5::DataSet data_aeq	     = distribution_file.getDataSet("aeq");
	h5::DataSet data_alpha 		 = distribution_file.getDataSet("alpha");
	h5::DataSet data_upar 		 = distribution_file.getDataSet("upar");
	h5::DataSet data_uper 		 = distribution_file.getDataSet("uper");
	h5::DataSet data_ppar 		 = distribution_file.getDataSet("ppar");
	h5::DataSet data_pper 		 = distribution_file.getDataSet("pper");
	h5::DataSet data_eta 		 = distribution_file.getDataSet("eta");
	h5::DataSet data_zeta 		 = distribution_file.getDataSet("zeta");
	h5::DataSet data_time 		 = distribution_file.getDataSet("time");
	h5::DataSet data_deta_dt     = distribution_file.getDataSet("deta_dt");
	h5::DataSet data_M_adiabatic = distribution_file.getDataSet("M_adiabatic");
	h5::DataSet data_Ekin 		 = distribution_file.getDataSet("Ekin");
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
	data_deta_dt.read(deta_dt_0);
	data_M_adiabatic.read(M_adiabatic_0);
	data_Ekin.read(Ekin_0);
	
	//Append to struct from single vector.
	for(int p=0; p<Constants::population; p++)
	{
		initial_particles[p].lamda.push_back(lamda_0.at(p));
		initial_particles[p].alpha.push_back(alpha_0.at(p));  
		initial_particles[p].aeq.push_back(aeq_0.at(p));
		initial_particles[p].ppar.push_back(ppar_0.at(p));
		initial_particles[p].pper.push_back(pper_0.at(p));
		initial_particles[p].upar.push_back(upar_0.at(p));
		initial_particles[p].uper.push_back(uper_0.at(p));
		initial_particles[p].Ekin.push_back(Ekin_0.at(p));
		initial_particles[p].time.push_back(time_0.at(p));
		initial_particles[p].zeta.push_back(zeta_0.at(p));
		initial_particles[p].eta.push_back(eta_0.at(p));
		initial_particles[p].deta_dt.push_back(deta_dt_0.at(p));
		initial_particles[p].M_adiabatic.push_back(M_adiabatic_0.at(p));
	}
	std::cout<<"\nParticle population: "<< initial_particles.size()<<std::endl;


//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
	int realthreads;   
	real wtime = omp_get_wtime();


	//---NOWPI---//
	if(argv[1]==string_no_wpi)
	{
		std::cout<<"\n\n"<<Constants::t_nowpi<<" sec NoWPI Simulation"<<std::endl;
		std::cout<<"\nExecution time estimation for 8 THREAD run: "<<(Constants::population*0.008/60) * Constants::t_nowpi <<" minutes."<<std::endl;
		std::cout<<"\nExecution time estimation for 20 THREAD run: "<<(Constants::population*0.017/60) * Constants::t_nowpi <<" minutes."<<std::endl;
		std::cout<<"\nForked..."<<std::endl;
		//---PARALLELISM Work sharing---//
		#pragma omp parallel
		{
			int id = omp_get_thread_num();
			if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
			#pragma omp for schedule(dynamic)
				for(int p=0; p<Constants::population; p++)     //dynamic because some chunks may have less workload.(particles can become invalid)
				{
					//Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
					no_wpi(p, lost_particles[p], ODPT, initial_particles[p]);
				}
		}	
		std::cout<<"\n"<<"Joined"<<std::endl;
		real time1 = omp_get_wtime()-wtime;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time1<<std::endl;
	}

	//Now initial_particles have the last state after NoWPI

	//---WPI---//
	if(argc==3)
	{
		//---LI---//
		if(argv[2]==string_li_wpi)
		{
			std::cout<<"\n\n"<<Constants::t_wpi<<" sec ray tracing WPI Simulation using Li formulas."<<std::endl;
			std::cout<<"Execution time estimation for 20 THREAD run: "<<(Constants::population*0.1/60) * Constants::t_wpi <<" minutes."<<std::endl;
			std::cout<<"\nForked..."<<std::endl;
			//---PARALLELISM Work sharing---//
			#pragma omp parallel
			{
				int id = omp_get_thread_num();
				if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
				#pragma omp for schedule(dynamic)
					for(int p=0; p<Constants::population; p++)     
					{
						//Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
						li_wpi(p, lost_particles[p], ODPT, initial_particles[p]);  //LI   + RAY TRACING
					}
			}	
			std::cout<<"\n"<<"Joined"<<std::endl;
			real time2 = omp_get_wtime()-wtime ;
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}

		//---BELL---//
		else if(argv[2]==string_bell_wpi)
		{
			std::cout<<"\n\n"<<Constants::t_wpi<<" sec WPI Simulation using Bell formulas. Wave magnitude(T): "<<Constants::By_wave<<std::endl;
			std::cout<<"Execution time estimation for 8 THREAD run: "<<(Constants::population*0.036/60) * Constants::t_wpi <<" minutes."<<std::endl;
			//---PARALLELISM Work sharing---//
			#pragma omp parallel
			{
				int id = omp_get_thread_num();
				if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
				#pragma omp for schedule(dynamic)
					for(int p=0; p<Constants::population; p++)     
					{
						//Void Function for particle's motion. Involves RK4 for Nsteps. Detected particles are saved in ODPT object, which is passed here by reference.
						bell_wpi(p, lost_particles[p], ODPT, initial_particles[p]); 		//BELL + THE WAVE IS EVERYWHERE
					}
			}	
			std::cout<<"\n"<<"Joined"<<std::endl;
			real time2 = omp_get_wtime()-wtime ;
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}
	}
//------------------------------------------------------------------ SIMULATION: END ---------------------------------------------------------------------//

//------------------------------------------------------------ OUTPUT DATA HDF5 --------------------------------------------------------------------------//
	//Assign from struct to vectors.
	std::vector<real> precip_id;
	std::vector<real> precip_lamda;
	std::vector<real> precip_alpha;
	std::vector<real> precip_aeq;
	std::vector<real> precip_time;
		
	for(int p=0; p<Constants::population; p++) 
	{
		if(!lost_particles[p].id.empty()) //Only precipitated were saved, other vectors in the struct are empty.
		{
			precip_id.push_back(lost_particles[p].id.at(0));
			precip_lamda.push_back(lost_particles[p].lamda.at(0)); 
			precip_alpha.push_back(lost_particles[p].alpha.at(0));
			precip_aeq.push_back(lost_particles[p].aeq.at(0));
			precip_time.push_back(lost_particles[p].time.at(0));
		}
	}

	h5::File file("h5files/test.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	//Detected particles
	h5::DataSet detected_lamda      = file.createDataSet("ODPT.lamda", ODPT.lamda);
	h5::DataSet detected_time       = file.createDataSet("ODPT.time", ODPT.time);
	h5::DataSet detected_id         = file.createDataSet("ODPT.id", ODPT.id);
	h5::DataSet detected_alpha      = file.createDataSet("ODPT.alpha", ODPT.alpha);
	h5::DataSet detected_aeq        = file.createDataSet("ODPT.aeq", ODPT.aeq);

	//Simulation data and Telescope specification - Scalars 
	h5::DataSet telescope_lamda    = file.createDataSet("ODPT.latitude", ODPT.latitude);
	h5::DataSet population         = file.createDataSet("population", 	Constants::population);
	h5::DataSet lamda_start_d      = file.createDataSet("lamda_start_d",Constants::lamda_start_d);
	h5::DataSet lamda_end_d        = file.createDataSet("lamda_end_d",  Constants::lamda_end_d);
	h5::DataSet aeq_start_d        = file.createDataSet("aeq_start_d",  Constants::aeq_start_d);
	h5::DataSet aeq_end_d          = file.createDataSet("aeq_end_d",    Constants::aeq_end_d);
	h5::DataSet Ekev0	           = file.createDataSet("Ekev0",   		Constants::Ekev0);
	h5::DataSet t			       = file.createDataSet("t", 			Constants::t);
	h5::DataSet By_wave            = file.createDataSet("By_wave",		Constants::By_wave);
	
	//Saved Particles that Precipitate.
	h5::DataSet saved_id     = file.createDataSet("precip_id", precip_id);
	h5::DataSet saved_lamda  = file.createDataSet("precip_lamda", precip_lamda);
	h5::DataSet saved_alpha  = file.createDataSet("precip_alpha", precip_alpha);
	h5::DataSet saved_aeq    = file.createDataSet("precip_aeq", precip_aeq);
	h5::DataSet saved_time   = file.createDataSet("precip_time", precip_time);


//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//
return 0; 

}


