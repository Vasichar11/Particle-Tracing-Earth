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
	std::cout<<"\n| "<<Constants::aeq_start_d << "  "<< " " << Constants::aeq_end_d <<"|\n";
	//std::cout<<"\nWith lamda distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\nWith latitude uniformly distributed, with varying interval"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::lamda_start_d << "  "<< " " << Constants::lamda_end_d <<"|\n";

	Particles single; //Single particle struct.
	std::vector<Particles> dstr(Constants::population, single);	//Vector of structs for particle distribution.
	
	real lamda0,Blam0;
	real lamda0_mr,Blam0_mr;
	real k,aeq0,salpha0,alpha0;
	real k_mr,aeq0_mr,salpha0_mr,alpha0_mr;
	real Beq0 = Bmag_dipole(0);   	 //Beq isn't always Beq0?

	std::random_device seed;         //Random seed. 
    std::mt19937 generator(seed());  //PRNG initialized with seed.
	real number;
	real interval;

	int p=0;
	int aeq_count = 0;
	while(aeq_count<Constants::aeq_dstr)
	{
		aeq0 = (Constants::aeq_start_d + aeq_count*Constants::aeq_step_d) * Constants::D2R;	  	  //linspace(start,stop,aeq_dstr)					
		//aeq0_mr = aeq0 + M_PI/2;		 													  	  //Mirror aeq  														
		//Varying interval for latitude selection - 2x2 linear system.
		//Aeq defines valid latitude interval.
		//Aeq~90, interval is small. Moving away from 90, interval increases.
		if      (0<aeq0 && aeq0<(M_PI/2)) 		interval = 180 - (2*aeq0*Constants::R2D); 
		else if ((M_PI/2)<aeq0 && aeq0<M_PI)    interval = (2*aeq0*Constants::R2D) - 180;
		//else if (aeq0==M_PI/2) 					interval = 0.1 ; //Equator particles, aeq=90, lamda~=0
		else if (aeq0==M_PI/2) 					interval = 0;    //Try only lamda=0?? Otherwise these particles would develop huge negative aeq for lamda=~0. 
		else { std::cout<<"\nError aeq0 initialization. aeq0= " << aeq0*Constants::R2D << std::endl; return EXIT_FAILURE;}

		int lat_count = 0;
		while(lat_count<Constants::lamda_dstr)												   //"Brute Force domain validation". Better solution? 
		{																				   	   //Loop until <lamda_dstr> valid states for this aeq0.
			std::uniform_real_distribution <real> distribution(-interval/2, interval/2);       //Uniform distribution with varying interval.
			number  = distribution(generator);											       //When aeq~90, interval is small. Moving away from 90, interval increases.
			lamda0 	= number * Constants::D2R;											    
			//lamda0_mr 	= - lamda0;											 			   //Mirror lamda

			//Find P.A at lamda0.
			Blam0 	   = Bmag_dipole(lamda0);
			//Blam0_mr   = Bmag_dipole(lamda0_mr);
			salpha0    = sin(aeq0)*sqrt(Blam0/Beq0); 										   //(2.20) Bortnik thesis
			//salpha0_mr    = sin(aeq0_mr)*sqrt(Blam0_mr/Beq0);  
			if(   !( (salpha0>1) || (salpha0<-1) || (salpha0==0) )   )  //|| (salpha0_mr>1) || (salpha0_mr<-1) || (salpha0_mr==0)
			{																			 	   //NOR of these should be true. Otherwise domain error.
				//Projecting aeq from alpha
				k       = ((aeq0*   Constants::R2D>90)    ? 1 : 0);     					   //kEN...(here k=0 or 1 ?)
				//k_mr    = ((aeq0_mr*Constants::R2D>90)    ? 1 : 0); 		
				alpha0  = pow(-1,k)*asin(salpha0)+k*M_PI;			 						   // sinx = a => x=(-1)^k * asin(a) + k*pi
				//alpha0_mr  = pow(-1,k_mr)*asin(salpha0_mr)+k_mr*M_PI;			 	
				dstr[p].initialize(Constants::eta0,aeq0,alpha0,lamda0,Constants::Ekev0,Blam0,0,0,0);
				//dstr[p+1].initialize(Constants::eta0,aeq0_mr,alpha0_mr,lamda0_mr,Constants::Ekev0,Blam0_mr,0,0,0);
				//Print initial state of particles.
				std::cout<<"\nParticle"<<p<<" aeq0: "<< aeq0*Constants::R2D <<", lamda0: "<< lamda0*Constants::R2D <<" gives alpha0: "<<alpha0*Constants::R2D<<std::endl;	
				//std::cout<<"\nParticle"<<p+1<<" aeq0: "<< aeq0_mr*Constants::R2D <<", lamda0: "<< lamda0_mr*Constants::R2D <<" gives alpha0: "<<alpha0_mr*Constants::R2D<<std::endl;
				p++;//p+=2;
				lat_count++;//lat_count+=2;
			}
		}	
		aeq_count++;
	}
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//

	//AEQ0 DISTRIBUTION CHECK
	const int sector_range = 15;
	const int view = 180;
	const int sectors = view/sector_range;
	std::array<int, sectors> aeq0_bins;
	std::fill(std::begin(aeq0_bins), std::end(aeq0_bins), 0); 		  //Initialize array elements with 0
	int sec;
	for(int p=0; p<Constants::population; p++)
	{
		sec = floor((dstr[p].aeq.at(0)*Constants::R2D)/sector_range); //Which sector has this particle?
		aeq0_bins.at(sec) ++; 										  //This sector has this particle
	}
	std::cout<<"\nEquatorial P.A Initialization: ";
	for(int sector=0;sector<sectors;sector++)					      //Print population of P.A bins
	{				
		std::cout<<"\naeq0 range: "<<sector*sector_range<< " - " <<(sector+1)*sector_range<< " has " << aeq0_bins.at(sector) << " particles.";
	}

//----------------------------------------WRITE TO HDF5 FILE------------------------------------//
	std::vector<real> lamda_dstr(Constants::population);
	std::vector<real> alpha_dstr(Constants::population);
	std::vector<real> aeq_dstr(Constants::population);
	std::vector<real> upar_dstr(Constants::population);
	std::vector<real> uper_dstr(Constants::population);
	std::vector<real> ppar_dstr(Constants::population);
	std::vector<real> pper_dstr(Constants::population);
	std::vector<real> eta_dstr(Constants::population);
	std::vector<real> M_adiabatic_dstr(Constants::population);
	std::vector<real> deta_dt_dstr(Constants::population);
	std::vector<real> time_dstr(Constants::population);
	std::vector<real> Ekin_dstr(Constants::population);
	std::vector<real> zeta_dstr(Constants::population);
	//Assign from struct to 1d vectors.
	for(int p=0; p<Constants::population; p++)
	{
		alpha_dstr[p]      = dstr[p].alpha.at(0);
		lamda_dstr[p]      = dstr[p].lamda.at(0);
		aeq_dstr[p]        = dstr[p].aeq.at(0);
		ppar_dstr[p]       = dstr[p].ppar.at(0);
		pper_dstr[p]       = dstr[p].pper.at(0);
		upar_dstr[p]       = dstr[p].upar.at(0);
		uper_dstr[p]       = dstr[p].uper.at(0);
		eta_dstr[p]        = dstr[p].eta.at(0);
		zeta_dstr[p]       = dstr[p].zeta.at(0);
		deta_dt_dstr[p]    = dstr[p].deta_dt.at(0);
		Ekin_dstr[p]       = dstr[p].Ekin.at(0);
		M_adiabatic_dstr[p]= dstr[p].M_adiabatic.at(0);
		time_dstr[p]       = dstr[p].time.at(0);
	}
	h5::File file("h5files/distribution.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

	h5::DataSet data_lat            = file.createDataSet("lat", lamda_dstr);
	h5::DataSet data_aeq            = file.createDataSet("aeq", aeq_dstr);
	h5::DataSet data_alpha          = file.createDataSet("alpha", alpha_dstr);
	h5::DataSet data_upar           = file.createDataSet("upar", upar_dstr);
	h5::DataSet data_uper           = file.createDataSet("uper", uper_dstr);
	h5::DataSet data_ppar           = file.createDataSet("ppar", ppar_dstr);
	h5::DataSet data_pper           = file.createDataSet("pper", pper_dstr);
	h5::DataSet data_eta            = file.createDataSet("eta", eta_dstr);
	h5::DataSet data_zeta           = file.createDataSet("zeta", zeta_dstr);
	h5::DataSet data_time           = file.createDataSet("time", time_dstr);
	h5::DataSet data_deta_dt        = file.createDataSet("deta_dt", deta_dt_dstr);
	h5::DataSet data_M_adiabatic    = file.createDataSet("M_adiabatic", M_adiabatic_dstr);
	h5::DataSet data_Ekin           = file.createDataSet("Ekin", Ekin_dstr);
	h5::DataSet aeq0bins            = file.createDataSet("aeq0_bins", aeq0_bins);

//----------------------------------------WRITE TO HDF5 FILE------------------------------------//

    return 0;
}
