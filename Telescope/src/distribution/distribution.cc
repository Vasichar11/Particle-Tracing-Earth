#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <random> 
#include <iomanip>  //For std::setprecision()

#include "headers/common.h"
#include "headers/struct_Particles.h"   
#include "headers/functions.h"
#include "headers/constants.h"

//HDF5 I/O
#include <highfive/H5File.hpp>          
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>


namespace h5 = HighFive;


int main()
{

//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//
	//Object for particles.
	std::cout<<"\n\nParticle population: " << Constants::population << "\n";
	//std::cout<<"\n\nEta distribution in degrees"<<"\n|From "<<" To|";
	//std::cout<<"\n| "<<Constants::eta_start_d << "  "<< " " << Constants::eta_end_d <<"|\n";
	std::cout<<"\nWith aeq fixed step distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::aeq_start_d << "  "<< " " << Constants::aeq_end_d <<"|";
	std::cout<<"\naeq step is: "<<Constants::aeq_step_d<< " degrees\n";

	Particles single; //Single particle struct.
	std::vector<Particles> dstr(Constants::population, single);	//Vector of structs for particle distribution.
	real lamda0,Blam0,aeq0,salpha0,alpha0,k,lamda_start_d,lamda_end_d,lamda_step_d;
	real Beq0 = Bmag_dipole(0);   	 //Beq isn't always Beq0?



	int p=0;
	int aeq_count = 0;
	while(aeq_count<Constants::aeq_dstr)
	{
		aeq0 = (Constants::aeq_start_d + aeq_count*Constants::aeq_step_d) * Constants::D2R;	  	  //linspace(start,stop,aeq_dstr)					
		lamda_end_d = 90; //default ending value
		//Find first valid lamda0 for this alpha0
		for(real lamda0_d=0;lamda0_d<90;lamda0_d+=Constants::h)
		{
			Blam0 	   = Bmag_dipole(lamda0_d*Constants::D2R);
			salpha0    = sin(aeq0)*sqrt(Blam0/Beq0); 
			if(   (salpha0>1) || (salpha0<-1) || (salpha0==0)   ) 
			{	//invalid lamda0_d
				std::cout<<"\nInvalid when lamda0_d = "<<lamda0_d<<" salpha0 = "<<salpha0;
				lamda_end_d = lamda0_d - Constants::h; //the last valid value
				break;
			}
		}
		//Now we have the first and last valid lamda0 for this aeq0
		lamda_start_d  = -lamda_end_d;			 				
		lamda_step_d   = (lamda_end_d - lamda_start_d)/(Constants::lamda_dstr-1); 	
		std::cout<<"\nFor aeq0 = "<<aeq0*Constants::R2D<<" degrees";
		std::cout<<"\nThe lamda fixed step distribution in degrees"<<"\n|From "<<" To|";
		std::cout<<"\n| "<<lamda_start_d << "  "<< " " << lamda_end_d <<"|";
		std::cout<<"\nlamda step is: "<<lamda_step_d<< " degrees\n";

		int lamda_count = 0;
		while(lamda_count<Constants::lamda_dstr)												   //"Brute Force domain validation". Better solution? 
		{																				   	   //Loop until <lamda_dstr> valid states for this aeq0.
			lamda0 = (lamda_start_d + lamda_count*lamda_step_d) * Constants::D2R;	  	  //linspace(start,stop,lamda_dstr)					
			Blam0 	   = Bmag_dipole(lamda0);
			salpha0    = sin(aeq0)*sqrt(Blam0/Beq0); 										   //(2.20) Bortnik thesis

			if(   !( (salpha0>1) || (salpha0<-1) || (salpha0==0) )   )  //|| (salpha0_mr>1) || (salpha0_mr<-1) || (salpha0_mr==0)
			{																			 	   //NOR of these should be true. Otherwise domain error.
				//Projecting aeq from alpha
				k       = ((aeq0*   Constants::R2D>90)    ? 1 : 0);     					   //kEN...(here k=0 or 1 ?)
				alpha0  = pow(-1,k)*asin(salpha0)+k*M_PI;			 						   // sinx = a => x=(-1)^k * asin(a) + k*pi
				dstr[p].initialize(Constants::eta0,aeq0,alpha0,lamda0,Constants::Ekev0,Blam0,0,0);

				//Print initial state of particles.
				std::cout<<"\nParticle"<<p<<" aeq0: "<< aeq0*Constants::R2D <<", lamda0: "<< lamda0*Constants::R2D <<" gives alpha0: "<<alpha0*Constants::R2D<<std::endl;	
				
				p++; lamda_count++;
			}
			else
			{
				std::cout<<"\nError.";
			}
		}	
		aeq_count++;
	}
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//

	//AEQ0 DISTRIBUTION CHECK
	const int sector_range = 5;
	const int view = 180;
	const int sectors = view/sector_range;
	std::array<int, sectors> aeq0_bins;
	std::fill(std::begin(aeq0_bins), std::end(aeq0_bins), 0); 		  //Initialize array elements with 0
	int sec;
	for(int p=0; p<Constants::population; p++)
	{
		sec = floor((dstr[p].aeq_init*Constants::R2D)/sector_range); //Which sector has this particle?
		aeq0_bins.at(sec) ++; 										  //This sector has this particle
	}
	std::cout<<"\nEquatorial P.A Initialization: ";
	for(int sector=0;sector<sectors;sector++)					      //Print population of P.A bins
	{				
		std::cout<<"\naeq0 range: "<<sector*sector_range<< " - " <<(sector+1)*sector_range<< " has " << aeq0_bins.at(sector) << " particles.";
	}

//----------------------------------------WRITE TO HDF5 FILE------------------------------------//
	
	std::vector<real> lamda_dstr(Constants::population), alpha_dstr(Constants::population), aeq_dstr(Constants::population), upar_dstr(Constants::population), uper_dstr(Constants::population), ppar_dstr(Constants::population), pper_dstr(Constants::population), eta_dstr(Constants::population), M_adiabatic_dstr(Constants::population), time_dstr(Constants::population), Ekin_dstr(Constants::population), zeta_dstr(Constants::population);

	//Assign from struct to 1d vectors.
	for(int p=0; p<Constants::population; p++)
	{
		alpha_dstr[p]      = dstr[p].alpha_init;
		lamda_dstr[p]      = dstr[p].lamda_init;
		aeq_dstr[p]        = dstr[p].aeq_init;
		ppar_dstr[p]       = dstr[p].ppar_init;
		pper_dstr[p]       = dstr[p].pper_init;
		upar_dstr[p]       = dstr[p].upar_init;
		uper_dstr[p]       = dstr[p].uper_init;
		eta_dstr[p]        = dstr[p].eta_init;
		zeta_dstr[p]       = dstr[p].zeta_init;
		Ekin_dstr[p]       = dstr[p].Ekin_init;
		M_adiabatic_dstr[p]= dstr[p].M_adiabatic_init;
		time_dstr[p]       = dstr[p].time_init;
	}
	h5::File file("h5files/test_new.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

	h5::DataSet data_lat            = file.createDataSet("lat", lamda_dstr);
	h5::DataSet data_aeq            = file.createDataSet("aeq", aeq_dstr);
	h5::DataSet data_alpha          = file.createDataSet("alpha", alpha_dstr);
	h5::DataSet data_upar           = file.createDataSet("upar", upar_dstr);
	h5::DataSet data_uper           = file.createDataSet("uper", uper_dstr);
	h5::DataSet data_ppar           = file.createDataSet("ppar", ppar_dstr);
	h5::DataSet data_pper           = file.createDataSet("pper", pper_dstr);
	h5::DataSet data_eta            = file.createDataSet("eta",  eta_dstr);
	h5::DataSet data_zeta           = file.createDataSet("zeta", zeta_dstr);
	h5::DataSet data_time           = file.createDataSet("time", time_dstr);
	h5::DataSet data_M_adiabatic    = file.createDataSet("M_adiabatic", M_adiabatic_dstr);
	h5::DataSet data_Ekin           = file.createDataSet("Ekin", Ekin_dstr);
	h5::DataSet aeq0bins            = file.createDataSet("aeq0_bins", aeq0_bins);

//----------------------------------------WRITE TO HDF5 FILE------------------------------------//

    return 0;
}
