//Standard library 
#include <cmath>                        
#include <iostream>                 
#include <vector>   
#include <cinttypes>   //int64_t type
#include <chrono>
#include <iomanip>    //setprecision()

//Function/Struct headers
#include "headers/wpi.h"
#include "headers/common.h"
#include "headers/constants.h"
#include "headers/struct_Particles.h"  
#include "headers/struct_Telescope.h"

//HDF5 I/O
#include <highfive/H5File.hpp>                                      
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>


int main()
{	
    //Position of the Particle Telescope.				
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell);	

//------------------------------------------------------ DISTRIBUTION OF PARTICLES--------------------------------------------------------//

	std::cout<<"\n\nParticle population: " << Constants::population;
	std::cout<<"\n\nEta distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::eta_start_d << "  "<< " " << Constants::eta_end_d <<"|\n";
	std::cout<<"\nWith aeq distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::aeq_start_d << "  "<< " " << Constants::aeq_end_d <<"|\n";
	std::cout<<"\nWith lamda distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::lamda_start_d << "  "<< " " << Constants::lamda_end_d <<"|\n";

	Particles single; //Single particle struct.
	std::vector<Particles> eql_dstr(Constants::population, single);	//Vector of structs for particle distribution.
						//equally distributed
	
	real eta0,aeq0,lamda0;
	real Beq0 = Bmag_dipole(0);	//Beq isn't always Beq0?


	for(int e=0, i=0; e<Constants::eta_dstr; e++)		
	{   																							 
		eta0 = (Constants::eta_start_d + e*Constants::eta_step_d) * Constants::D2R;			   														
		
		for(int a=0; a<Constants::aeq_dstr; a++)											
		{	
			aeq0 = (Constants::aeq_start_d + a*Constants::aeq_step_d) * Constants::D2R;	         	  	    

			for(int l=0; l<Constants::lamda_dstr; l++, i++)
		    {
		    	lamda0 = (Constants::lamda_start_d + l*Constants::lamda_step_d) * Constants::D2R; 
			
				//Now push-back initial values in RADS.
			    eql_dstr[i].set_eta(eta0);  		
			    eql_dstr[i].set_aeq(aeq0);  
			    eql_dstr[i].set_lamda(lamda0);  
				
				//And push-back other resulting values.(P.A, speed, momentum)
				//Exceptions are involved
				try
				{	
					eql_dstr[i].calculations0(Beq0,lamda0,0,0,aeq0,Constants::Ekev0);
				}					
				catch(int exception)	 
				{
				
					if(exception == 99)
					{
						std::cout<< "\n" << "Caught exception: |sinusoidal|>1. P.A out of domain.\n" ;
						return EXIT_FAILURE;	
					}
					else if(exception == 98)	
					{
						std::cout<< "\n" << "Caught exception: aeq0=0||180.\n" ;
						return EXIT_FAILURE;
					}	
				
				}		
				//Print initial state of particles.
				std::cout<<"\nParticle"<<i<<" with eta0: "<< eta0*Constants::R2D <<", aeq0: "<< aeq0*Constants::R2D <<" and lamda0: "<<lamda0*Constants::R2D<< " in degrees.\n";				
			}
		}	
	}

//------------------------------------------------------ DISTRIBUTION OF PARTICLES: END---------------------------------------------------//
    
    std::cout.precision(8);    //Output precision of 16 decimals
    std::cout<<std::scientific; //Representation with e notation

//---------------------------------------------------------------SIMULATION---------------------------------------------------------------//
	auto rk_start = std::chrono::high_resolution_clock::now();
	for(int p=0; p<Constants::population; p++)     //Loop for all particles
	{
		//Void Function for particle's motion. Involves RK4 for Nsteps. 
        wpi(p, eql_dstr[p].lamda.front(), eql_dstr[p].alpha.front(), eql_dstr[p].aeq.front(), eql_dstr[p].ppar.front(), eql_dstr[p].pper.front(), eql_dstr[p].upar.front(), eql_dstr[p].uper.front(), eql_dstr[p].zeta.front(), eql_dstr[p].M_adiabatic.front(), eql_dstr[p].eta.front(), eql_dstr[p].time.front(), ODPT);
		//Detected particles are saved in ODPT object, which is passed here by reference.
	}	
	auto rk_stop = std::chrono::high_resolution_clock::now();  
	auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(rk_stop - rk_start);
	real rk_time = duration2.count() * pow(10,-6);
	std::cout << "\nExecution time " << std::setprecision(8) << rk_time <<" seconds\n" ;
//--------------------------------------------------------------SIMULATION: END------------------------------------------------------------//



//---------------------------------------------------------- OUTPUT DATA HDF5 -------------------------------------------------------------//


//--------------------------------------------FOR DETECTED PARTICLES------------------------------------------------//	
	//Create hdf5 file for particles that cross the satellite.
	h5::File file2("h5files/detected.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

	h5::DataSet dataset2_lamda = file2.createDataSet("detected_lamda 1D", ODPT.lamda);
	h5::DataSet dataset2_time  = file2.createDataSet("detected_time 1D",  ODPT.time);
	h5::DataSet dataset2_id    = file2.createDataSet("detected_id 1D",  ODPT.id);
	h5::DataSet dataset2_aeq   = file2.createDataSet("detected_aeq 1D", ODPT.aeq);
	h5::DataSet dataset2_alpha = file2.createDataSet("detected_alpha 1D", ODPT.alpha);
	h5::DataSet telescope_lamda= file2.createDataSet("telescope_lamda scalar", ODPT.latitude);
	h5::DataSet population     = file2.createDataSet("whole population", Constants::population);
	h5::DataSet initial_lamda  = file2.createDataSet("initial latitude", Constants::lamda0);
	h5::DataSet dataset_time   = file2.createDataSet("simulation time", Constants::t);

//----------------------------------------FOR ALL PARTICLES (IF NEEDED) --------------------------------------------//	
/*
	std::vector<std::vector<real>> aeq2_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>  aeq_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>  eta_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>lamda_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> time_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> upar_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	//ASSIGN FROM STRUCT INTO 2D VECTORS. BETTER WAY?
	for(int p=0; p<Constants::population; p++)
	{
	    for(unsigned long int i=0; i<eql_dstr[p].lamda.size(); i++)  //Fill 2D vector from vector of structs.
	    {
	        aeq_plot[p][i]  = eql_dstr[p].aeq.at(i);  
	        eta_plot[p][i]  = eql_dstr[p].eta.at(i);  
	        lamda_plot[p][i]= eql_dstr[p].lamda.at(i);  
	    	time_plot[p][i] = eql_dstr[p].time.at(i);
	    	upar_plot[p][i] = eql_dstr[p].upar.at(i);
	    }
	}
	//Create hdf5 file.
	h5::File file1("particles.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	//Create datasets and write vectors to them(shortcut syntax).
	h5::DataSet dataset1_lamda =  file1.createDataSet("lamda 2D",       	lamda_plot);
	h5::DataSet dataset1_aeq   =  file1.createDataSet("aeq 2D",         	aeq_plot);
	h5::DataSet dataset1_eta   =  file1.createDataSet("eta 2D",         	eta_plot);
	h5::DataSet dataset1_time  =  file1.createDataSet("simulation time",	time_plot);
	h5::DataSet dataset1_upar  =  file1.createDataSet("upar",	    		upar_plot);
*/
//----------------------------------------------------- OUTPUT DATA HDF5 : END ----------------------------------------------------------//


    std::cout << "\n\n";
    return 0; 

}




