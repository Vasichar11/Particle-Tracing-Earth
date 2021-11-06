#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <iomanip>  //For std::setprecision()
#include <omp.h>

//Same directory headers							    
//Preprocessor macro instructions are added in files to obey ODR.
#include "headers/adiabatic_motion.h"
#include "headers/common.h"
#include "headers/struct_Particles.h"   		    	
#include "headers/struct_Telescope.h"  				
#include "headers/constants.h"

//HDF5 I/O
#include <highfive/H5File.hpp>          
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;


int main()
{
	//Position of the Particle Telescope.		
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell);		
	
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//
	//Object for particles.
	std::cout<<"\n\nParticle testing population: " << Constants::test_pop << "\n";
	std::cout<<"\n\nEta distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::eta_start_d << "  "<< " " << Constants::eta_end_d <<"|\n";
	std::cout<<"\nWith aeq distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::aeq_start_d << "  "<< " " << Constants::aeq_end_d <<"|\n";
	std::cout<<"\nWith lamda distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::lamda_start_d << "  "<< " " << Constants::lamda_end_d <<"|\n";

	Particles single; //Single particle struct.
	std::vector<Particles> eql_dstr(Constants::test_pop, single);	//Vector of structs for particle distribution.
						//equally distributed
	
	real eta0,aeq0,lamda0,k;
	real Blam0,salpha0,alpha0;
	real Beq0 = Bmag_dipole(0);	//Beq isn't always Beq0?

	for(int e=0, p=0; e<Constants::eta_dstr; e++)		
	{   																							 
		eta0 = (Constants::eta_start_d + e*Constants::eta_step_d) * Constants::D2R;			   														
		
		for(int a=0; a<Constants::aeq_dstr; a++)											
		{	
			aeq0 = (Constants::aeq_start_d + a*Constants::aeq_step_d) * Constants::D2R;	         	  	    

			for(int l=0; l<Constants::lamda_dstr; l++, p++)
		    {
		    	lamda0 = (Constants::lamda_start_d + l*Constants::lamda_step_d) * Constants::D2R;  	
				
				//Find P.A at lamda0.
				Blam0=Bmag_dipole(lamda0);
				salpha0=sin(aeq0)*sqrt(Blam0/Beq0); //(2.20) Bortnik thesis
				if( (salpha0<-1) || (salpha0>1) || (salpha0==0) ) {eql_dstr.pop_back(); p--; continue; } //Exclude these particles.
				k = ((aeq0*Constants::R2D>90) ? 1 : 0);
				alpha0=pow(-1,k)*asin(salpha0)+k*M_PI;		//If aeq0=150 => alpha0=arcsin(sin(150))=30 for particle in equator.Distribute in alpha instead of aeq?		 	
				eql_dstr[p].initialize(eta0,aeq0,alpha0,lamda0,Constants::Ekev0,Blam0,0,0,0);

				//Print initial state of particles.
				//std::cout<<"\nParticle"<<p<<" aeq0: "<< aeq0*Constants::R2D <<", lamda0: "<< lamda0*Constants::R2D <<" gives alpha0: "<<alpha0*Constants::R2D<<std::endl;				
			}
		}	
	}
	int64_t track_pop = eql_dstr.size(); //Population of particles that will be tracked.
	std::cout<<"\n"<<Constants::test_pop - track_pop<<" Particles were excluded from the initial population due to domain issues.\nThe particle population for the tracer is now: "<<track_pop<<"\n";
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//
	
    
//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
    int realthreads;   
	//---PARALLELISM Work sharing---//
	double wtime = omp_get_wtime();
	std::cout<<"\nForked..."<<std::endl;
 	#pragma omp parallel
    {
    	int id = omp_get_thread_num();
		if(id==0){ realthreads = omp_get_num_threads();}
		
		#pragma omp for schedule(static)
			for(int p=0; p<track_pop; p++)     //Chunk=1, pass blocks of 1 iteration to each thread.
			{
				//std::cout<<"\nBouncing particle "<<p<<" "<<id<<std::flush;
				//Void Function for particle's motion. Involves RK4 for Nsteps. 
				//Detected particles are saved in ODPT object, which is passed here by reference.
				adiabatic_motion(p, eql_dstr[p], ODPT);   
			}	
	}
    std::cout<<"\n"<<"Joined"<<std::endl;
	wtime = omp_get_wtime()-wtime;
	std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<wtime<<std::endl;
//------------------------------------------------------------------ SIMULATION: END ---------------------------------------------------------------------//



//------------------------------------------------------------ OUTPUT DATA HDF5 --------------------------------------------------------------------------//

	std::vector<std::vector<real>>  aeq_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> lamda_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	
	//std::vector<std::vector<real>> time_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	//std::vector<std::vector<real>> alpha_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	//std::vector<std::vector<real>> deta_dt_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	

	//Assign from struct to 2d vectors.
	for(int p=0; p<track_pop; p++)
	{
	    //for(size_t i=0; i<eql_dstr[p].alpha.size(); i++)  
	    //{
	    	//time_plot[p][0] = eql_dstr[p].time.at(0);  //take the initial values
	    	aeq_plot[p][0] = eql_dstr[p].aeq.at(0);
	    	//alpha_plot[p][i] = eql_dstr[p].alpha.at(i);
	    	lamda_plot[p][0] = eql_dstr[p].lamda.at(0);
	    	//deta_dt_plot[p][i] = eql_dstr[p].deta_dt.at(i);
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
	
	//Saved Particles
	h5::DataSet saved_lamda  = file.createDataSet("lamda_plot", lamda_plot);
	//h5::DataSet saved_alpha  = file.createDataSet("alpha_plot", alpha_plot);
	//h5::DataSet saved_deta   = file.createDataSet("deta_dt", deta_dt_plot);
	h5::DataSet saved_aeq    = file.createDataSet("aeq_plot", aeq_plot);
	//h5::DataSet saved_time   = file.createDataSet("time_plot", time_plot);




//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//

return 0; 

}


