#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <array> 
#include <random> 
#include <iomanip>  //For std::setprecision()
#include <omp.h>

//Same directory headers							    
//Preprocessor macro instructions are added in files to obey ODR.
#include "headers/bell_nowpi.h"
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
	if(argc!=2)
	{
		std::cout<<"Error. Argc should be 2. Set second argv from the list:(nowpi, wpi, wpi_ray)."<<std::endl;
		return EXIT_FAILURE;
	}

	
	//Position of the Particle Telescope.		
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell);		
	
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//
	//Object for particles.
	std::cout<<"\n\nParticle testing population: " << Constants::test_pop << "\n";
	//std::cout<<"\n\nEta distribution in degrees"<<"\n|From "<<" To|";
	//std::cout<<"\n| "<<Constants::eta_start_d << "  "<< " " << Constants::eta_end_d <<"|\n";
	std::cout<<"\nWith aeq fixed step distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::aeq_start_d << "  "<< " " << Constants::aeq_end_d <<"|\n";
	//std::cout<<"\nWith lamda distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\nWith latitude uniformly distributed, with varying interval"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::lamda_start_d << "  "<< " " << Constants::lamda_end_d <<"|\n";

	Particles single; //Single particle struct.
	std::vector<Particles> dstr(Constants::test_pop, single);	//Vector of structs for particle distribution.
	
	real lamda0,Blam0;
	real k,aeq0,salpha0,alpha0;
	real Beq0 = Bmag_dipole(0);	//Beq isn't always Beq0?

	std::random_device seed;         //random seed 
    std::mt19937 generator(seed());  //PRNG initialized with seed
	real number, count;
	real interval;


	for(int a=0, p=0; a<Constants::aeq_dstr; a++)
	{
		aeq0 = (Constants::aeq_start_d + a*Constants::aeq_step_d) * Constants::D2R;	//Distributed with fixed step in fixed interval.		   														
		
		//Varying interval for latitude selection - 2x2 linear system.
		
		if (0<aeq0 && aeq0<90) 		  interval = 180 - aeq0;
		
		else if (90<aeq0 && aeq0<180) interval = aeq0 - 180;
		
		else { std::cout<<"\nError aeq0 initialization."; return EXIT_FAILURE;}

		//"Brute Force domain validation". Better solution? 
		//Loop until <lamda_dstr> valid states for this aeq0 -> Flat aeq distribution * Weighting factors = Realistic distribution
		count = 0;
		while(count<Constants::lamda_dstr) 
		{
			std::uniform_real_distribution <real> distribution(-interval/2, interval/2); //Uniform distribution with varying interval.
			number  = distribution(generator);											 //lat~90, interval is small. Moving away from 90, interval increases.
			lamda0 	= number * Constants::D2R; 											 //Done to increase particle randomness in latitude. 
			//Find P.A at lamda0.
			Blam0 	   = Bmag_dipole(lamda0);
			salpha0    = sin(aeq0)*sqrt(Blam0/Beq0);  //(2.20) Bortnik thesis
			if(   !((salpha0>1) || (salpha0<-1) || (salpha0==0))   ) //NOR of these should be true. Otherwise domain error. 
			{
				//Projecting aeq from alpha
				k    = ((aeq0*Constants::R2D>90)    ? 1 : 0); 		//kEN...(here k=0 or 1 ?)
				alpha0=pow(-1,k)*asin(salpha0)+k*M_PI;			 	// sinx = a => x=(-1)^k * asin(a) + k*pi
				dstr[p].initialize(Constants::eta0,aeq0,alpha0,lamda0,Constants::Ekev0,Blam0,0,0,0);
				p ++ ;
				count++;
				//Print initial state of particles.
				//std::cout<<"\nParticle"<<p<<" aeq0: "<< aeq0*Constants::R2D <<", lamda0: "<< lamda0*Constants::R2D <<" gives alpha0: "<<alpha0*Constants::R2D<<std::endl;	
			}
		}	
	}

	int64_t track_pop = dstr.size(); //Population of particles that will be tracked.
	std::cout<<"\n"<<Constants::test_pop - track_pop<<" Particles were excluded from the initial population due to domain issues.\nThe particle population for the tracer is now: "<<track_pop<<"\n";
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//
	
	//AEQ0 DISTRIBUTION
	const int sector_range = 15;
	const int view = 180;
	const int sectors = view/sector_range;
	std::array<int, sectors> aeq0_bins;
	int sec;
	for(int p=0; p<track_pop; p++)
	{
		sec = floor((dstr[p].aeq.at(0)*Constants::R2D)/sector_range);
		aeq0_bins.at(sec) ++;
	}
	for(int sector=0;sector<sectors;sector++)
	{
		std::cout<<"\nSector: "<<sector<< " has " << aeq0_bins.at(sector) << " particles.";
	}
//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
	int realthreads;   
	//---PARALLELISM Work sharing---//
	real wtime = omp_get_wtime();
	std::string s1("nowpi");
	std::string s2("wpi");
	std::string s3("wpi_ray");

	if( !(s1.compare(argv[1])) )
	{
        std::cout<<"\n\nNoWPI Simulation using Bell formulas"<<std::endl;
		std::cout<<"\nForked..."<<std::endl;
		#pragma omp parallel
    	{

    	int id = omp_get_thread_num();
		if(id==0){ realthreads = omp_get_num_threads();}
		
		#pragma omp for schedule(dynamic)
			for(int p=0; p<track_pop; p++)     //dynamic because some chunks may have less workload.(particles can become invalid)
			{
				//std::cout<<"\nBouncing particle "<<p<<" "<<id<<std::flush;
				//Void Function for particle's motion. Involves RK4 for Nsteps. 
				//Detected particles are saved in ODPT object, which is passed here by reference.
				nowpi(p, dstr[p], ODPT);
			}
		}	
    	std::cout<<"\n"<<"Joined"<<std::endl;
		wtime = omp_get_wtime()-wtime;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<wtime<<std::endl;
	}
	
	else if( !(s2.compare(argv[1])) )
	{
        std::cout<<"\n\nWPI Simulation using Bell formulas. Wave magnitude(T): "<<Constants::By_wave<<std::endl;
		std::cout<<"\nForked..."<<std::endl;
		#pragma omp parallel
    	{

    	int id = omp_get_thread_num();
		if(id==0){ realthreads = omp_get_num_threads();}
		
		#pragma omp for schedule(dynamic)
			for(int p=0; p<track_pop; p++)     //dynamic because some chunks may have less workload.(particles can become invalid)
			{
				//std::cout<<"\nBouncing particle "<<p<<" "<<id<<std::flush;
				//Void Function for particle's motion. Involves RK4 for Nsteps. 
				//Detected particles are saved in ODPT object, which is passed here by reference.
				wpi(p, dstr[p], ODPT);
			}
		}	
    	std::cout<<"\n"<<"Joined"<<std::endl;
		wtime = omp_get_wtime()-wtime;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<wtime<<std::endl;
	}
	
	else if( !(s3.compare(argv[1])) )
	{
        std::cout<<"\n\nWPI Simulation using Li formulas. Read Ray from CSV"<<std::endl;
		std::cout<<"\nForked..."<<std::endl;
		#pragma omp parallel
    	{

    	int id = omp_get_thread_num();
		if(id==0){ realthreads = omp_get_num_threads();}
		
		#pragma omp for schedule(dynamic)
			for(int p=0; p<track_pop; p++)     //dynamic because some chunks may have less workload.(particles can become invalid)
			{
				//std::cout<<"\nBouncing particle "<<p<<" "<<id<<std::flush;
				//Void Function for particle's motion. Involves RK4 for Nsteps. 
				//Detected particles are saved in ODPT object, which is passed here by reference.
				wpi_ray(p, dstr[p], ODPT);
			}
		}	
    	std::cout<<"\n"<<"Joined"<<std::endl;
		wtime = omp_get_wtime()-wtime;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<wtime<<std::endl;
	}

	else
	{
		std::cout<<"\nArgument variable doesn't match any of the program's possible implementations.\n\nTry nowpi, wpi, or wpi_ray as the second argument variable."<<std::endl;
		return EXIT_FAILURE;	
	}
//------------------------------------------------------------------ SIMULATION: END ---------------------------------------------------------------------//


//------------------------------------------------------------ OUTPUT DATA HDF5 --------------------------------------------------------------------------//

	std::vector<real>  aeq0_plot(track_pop);
	std::vector<real> lamda0_plot(track_pop);
	
	//std::vector<std::vector<real>> time_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	//std::vector<std::vector<real>> alpha_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	//std::vector<std::vector<real>> deta_dt_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	

	//Assign from struct to 2d vectors.
	for(int p=0; p<track_pop; p++)
	{
	    //for(size_t i=0; i<dstr[p].alpha.size(); i++)  
	    //{
	    	//time_plot[p][0] = dstr[p].time.at(0);  //take the initial values
	    	aeq0_plot[p] = dstr[p].aeq.at(0);
	    	//alpha_plot[p][i] = dstr[p].alpha.at(i);
	    	lamda0_plot[p] = dstr[p].lamda.at(0);
	    	//deta_dt_plot[p][i] = dstr[p].deta_dt.at(i);
	    //}
	}
    
	h5::File file("h5files/detected.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	//Detected particles
	h5::DataSet detected_lamda      = file.createDataSet("ODPT.lamda", ODPT.lamda);
	h5::DataSet detected_time       = file.createDataSet("ODPT.time", ODPT.time);
	h5::DataSet detected_id         = file.createDataSet("ODPT.id", ODPT.id);
	h5::DataSet detected_alpha      = file.createDataSet("ODPT.alpha", ODPT.alpha);

	//Simulation data and Telescope specification - Scalars 
	h5::DataSet telescope_lamda    = file.createDataSet("ODPT.latitude", ODPT.latitude);
	h5::DataSet population         = file.createDataSet("population", track_pop);
	h5::DataSet lamda_start_d      = file.createDataSet("lamda_start_d",Constants::lamda_start_d);
	h5::DataSet lamda_end_d        = file.createDataSet("lamda_end_d",  Constants::lamda_end_d);
	h5::DataSet aeq_start_d        = file.createDataSet("aeq_start_d",  Constants::aeq_start_d);
	h5::DataSet aeq_end_d          = file.createDataSet("aeq_end_d",    Constants::aeq_end_d);
	h5::DataSet Ekev0	           = file.createDataSet("Ekev0",   		Constants::Ekev0);
	h5::DataSet t			       = file.createDataSet("t", Constants::t);
	h5::DataSet By_wave            = file.createDataSet("By_wave",Constants::By_wave);
	h5::DataSet aeq0bins           = file.createDataSet("aeq0_bins", aeq0_bins);
	
	//Saved Particles
	h5::DataSet saved_lamda0  = file.createDataSet("lamda0_plot", lamda0_plot);
	//h5::DataSet saved_alpha  = file.createDataSet("alpha_plot", alpha_plot);
	//h5::DataSet saved_deta   = file.createDataSet("deta_dt", deta_dt_plot);
	h5::DataSet saved_aeq0    = file.createDataSet("aeq0_plot", aeq0_plot);
	//h5::DataSet saved_time   = file.createDataSet("time_plot", time_plot);


//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//

return 0; 

}


