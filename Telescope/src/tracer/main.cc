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
	if( (argc<2 || argc>3)   ||    argv[1]!=string_no_wpi    ||    ( argc==3 && argv[2]!=string_bell_wpi && argv[2]!=string_li_wpi ) )
	{ 
		std::cout<<"\nArgument variables don't match any of the program's possible implementations.\nSet second argv=no_wpi and third argv from the list:(bell_wpi, li_wpi) to introduce a wave.\n"<<std::endl;
		return EXIT_FAILURE;	
	}


//------------------------------------------------------------READ AND ASSIGN DISTRIBUTION FROM H5 FILE --------------------------------------------------------------//
	h5::File distribution_file("h5files/1p_testAEQ_testLAMDA.h5", h5::File::ReadOnly);
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
	if(argv[1]==string_no_wpi)
	{
		std::cout<<"\n\n"<<Constants::t_nowpi<<" sec NoWPI Simulation"<<std::endl;
		std::cout<<"\nExecution time estimation for 8 THREAD run: "<<(Population*0.008/60) * Constants::t_nowpi <<" minutes."<<std::endl;
		std::cout<<"Execution time estimation for 20 THREAD run: "<<(Population*0.017/60) * Constants::t_nowpi <<" minutes."<<std::endl;
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
					no_wpi(p, dstr[p], ODPT);
				}
		}	
		std::cout<<"\n\n"<<"Joined"<<std::endl;
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
			std::cout<<"Execution time estimation for 20 THREAD run: "<<(Population*0.1/60) * Constants::t_wpi <<" minutes."<<std::endl;
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
						li_wpi(p, dstr[p], ODPT);  //LI   + RAY TRACING
					}
			}	
			std::cout<<"\n\n"<<"Joined"<<std::endl;
			real time2 = omp_get_wtime()-wtime ;
			std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
		}

		//---BELL---//
		else if(argv[2]==string_bell_wpi)
		{
			std::cout<<"\n\n"<<Constants::t_wpi<<" sec WPI Simulation using Bell formulas. Wave magnitude(T): "<<Constants::By_wave<<std::endl;
			std::cout<<"Execution time estimation for 8 THREAD run: "<<(Population*0.036/60) * Constants::t_wpi <<" minutes."<<std::endl;
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
						bell_wpi(p, dstr[p], ODPT); 		//BELL + THE WAVE IS EVERYWHERE
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
	std::vector<real> precip_id, precip_lamda, precip_alpha, precip_aeq, precip_time, neg_id, nan_id, high_id, lamda00, ppar00, pper00, alpha00, aeq00, eta00, time00;

	for(int p=0; p<Population; p++) 
	{
		//Last particle states(that can become first states for next simulation).
		lamda00.push_back(dstr[p].lamda00);
    	ppar00.push_back(dstr[p].ppar00); 
    	pper00.push_back(dstr[p].pper00); 
    	alpha00.push_back(dstr[p].alpha00); 
    	aeq00.push_back(dstr[p].aeq00); 
    	eta00.push_back(dstr[p].eta00); 
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
	}




	//File name based on the argument variables
	std::string  save_file = "h5files/" + std::to_string(Constants::population) + "p_";   
	if 	  	(argc==2) save_file = save_file  + std::string(argv[1]) ;
	else if (argc==3) save_file = save_file  + std::string("both"); 
	save_file = save_file + ".h5";

	h5::File file(save_file, h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	//Simulation data and Telescope specification - Scalars 
	h5::DataSet telescope_lamda= file.createDataSet("ODPT.latitude", ODPT.latitude);
	h5::DataSet data_population= file.createDataSet("population", 	Population);
	h5::DataSet Ekev0	       = file.createDataSet("Ekev0",   		Constants::Ekev0);
	h5::DataSet t			   = file.createDataSet("t", 			Constants::t);
	
	//Detected particles
	h5::DataSet detected_lamda = file.createDataSet("ODPT.lamda", ODPT.lamda);
	h5::DataSet detected_time  = file.createDataSet("ODPT.time", ODPT.time);
	h5::DataSet detected_id    = file.createDataSet("ODPT.id", ODPT.id);
	h5::DataSet detected_alpha = file.createDataSet("ODPT.alpha", ODPT.alpha);
	h5::DataSet detected_aeq   = file.createDataSet("ODPT.aeq", ODPT.aeq);

	//Precipitating Particles
	h5::DataSet precipitated_id        = file.createDataSet("precip_id", 	 precip_id);
	h5::DataSet precipitated_lamda     = file.createDataSet("precip_lamda",  precip_lamda);
	h5::DataSet precipitated_aeq       = file.createDataSet("precip_aeq", 	 precip_aeq);
	h5::DataSet precipitated_time      = file.createDataSet("precip_time",	 precip_time);

	//Negative P.A Particles
	h5::DataSet data_neg_id           = file.createDataSet("neg_id", 	 neg_id);

	//High P.A Particles
	h5::DataSet data_high_id          = file.createDataSet("high_id", 	 high_id);

	//High P.A Particles
	h5::DataSet data_nan_id           = file.createDataSet("nan_id", 	 nan_id);


	//Particles states at noWPI end.
	h5::DataSet ending_lamda   = file.createDataSet("lamda00", lamda00);
	h5::DataSet ending_ppar    = file.createDataSet("ppar00",  ppar00);
	h5::DataSet ending_pper    = file.createDataSet("pper00",  pper00);
	h5::DataSet ending_alpha   = file.createDataSet("alpha00", alpha00);
	h5::DataSet ending_aeq     = file.createDataSet("aeq00",   aeq00);
	h5::DataSet ending_eta     = file.createDataSet("eta00",   eta00);
	h5::DataSet ending_time    = file.createDataSet("time00",  time00);



//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//
return 0; 

}


