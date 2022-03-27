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




//Member function to find valid latitude range
void lamda_domain(real aeq0, real &lamda_start_d, real &lamda_end_d)
{
	const real Beq0 = Bmag_dipole(0);   	 //Beq isn't always Beq0?
	lamda_end_d     = 0;
	lamda_start_d   = 0;
	real salpha0;
	do
	{
		//("Brute Force domain validation")
		real Blam0   = Bmag_dipole(lamda_end_d*Constants::D2R);
		salpha0      = sin(aeq0)*sqrt(Blam0/Beq0);  //salpha = sin(aeq)*sqrt(Blam/Beq)
		lamda_end_d += Constants::h; //Gradually increase and check if salpha0 is valid.
	}
	while( (salpha0<=1) && (salpha0>=-1) && (salpha0!=0) ) ; // &&lamda0<M_PI
	lamda_end_d   -= Constants::h;     //This is the last (positive) valid value.
	lamda_start_d  = -lamda_end_d;	 //This would be the first (negative) valid value. 		
}



int main(int argc, char **argv)
{

	//std::cout.precision(8);			//Output 16 decimal precise
	//std::cout<<std::scientific;		//For e notation representation
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//
	
	std::string string_evenly   = "evenly";
	std::string string_normal   = "normal";
	std::string string_uniform  = "uniform";
	std::string string_test     = "test";
	
	//---ARGV ERROR---//
	if( argc!=3   ||   ( argv[1]!=string_evenly  &&  argv[1]!=string_uniform  &&  argv[1]!=string_normal  &&   argv[1]!=string_test)   ||   ( argv[2]!=string_evenly  &&   argv[2]!=string_uniform   &&  argv[2]!=string_normal  &&   argv[2]!=string_test) )
	{ 
		std::cout<<"\nArgument variables don't match any of the program's possible implementations.\nSet second and third argvs from the list:(evenly, uniform, normal, test) to distribute the particles.\n"<<std::endl;
		return EXIT_FAILURE;	
	}
	std::cout<<"\n\nParticle population: " << Constants::population;
	//std::cout<<"\n\nEta distribution in degrees"<<"\n|From "<<" To|";
	//std::cout<<"\n| "<<Constants::eta_start_d << "  "<< " " << Constants::eta_end_d <<"|\n";


	Particles single; //Single particle struct.
	std::vector<Particles> dstr(Constants::population, single);	//Vector of structs for particle distribution.
	real lamda0,Blam0,aeq0,salpha0,alpha0,k;
	//const real Beq0 = Bmag_dipole(0);   	 //Beq isn't always Beq0?
	std::random_device seed;         //Random seed. 
	std::mt19937 generator(seed());  //PRNG initialized with seed.
	int p=0;
	int aeq_count = 0;
	//Loop for <lamda_dstr> different particle latitudes.
	while(aeq_count<Constants::aeq_dstr)  
	{
		//----------------P.A dstr-----------------//
		//Evenly dstr
		if(argv[1]==string_evenly)   
		{
			aeq0 = (Constants::aeq_start_d + aeq_count*Constants::aeq_step_d) * Constants::D2R;	  	  //linspace(start,stop,aeq_dstr). Goal is: symmetry				
		}
		//Uniformly dstr
		else if (argv[1]==string_uniform)   			
		{
			std::uniform_real_distribution <real> uniform_aeq  (Constants::aeq_start_d, Constants::aeq_end_d);     //All values in this range are equally probable. Goal is: randomness
			aeq0  = uniform_aeq(generator) * Constants::D2R; 	
		}
		//Normaly dstr
		else if (argv[1]==string_normal)   			
		{
			real number;
			do
			{
				std::normal_distribution <real> normal_aeq(Constants::mean_aeq, Constants::stdev_aeq);		   //Goal is: reach normal dstr in logarithm scale.
				number = normal_aeq(generator);

			}while(number<=0 || number>=180); //Until valid aeq=(0,180).
			aeq0  = number * Constants::D2R; 	
		}
		else if (argv[1]==string_test)
		{
			aeq0 = Constants::aeq0;
		}
		else return EXIT_FAILURE;
		//----------------P.A dstr-----------------//


		//-------------Latitude domain-------------//
		real lamda_end_d, lamda_start_d;
		lamda_domain(aeq0, lamda_start_d, lamda_end_d); //Finds lamda_start_d && lamda_end_d. 
		//std::cout<<"\nFor aeq0 = "<<aeq0*Constants::R2D<<" degrees\nThe latitude domain in degrees"<<"\n|From "<<" To|\n| "<<lamda_start_d << "  "<< " " << lamda_end_d <<"|\n";
		//-------------Latitude domain-------------//


		int lamda_count = 0;
		while(lamda_count<Constants::lamda_dstr)  
		{
		//--------------Latitude dstr--------------//
			if(argv[2]==string_evenly)   
			{
				lamda0 = (Constants::lamda_start_d + lamda_count*Constants::lamda_step_d) * Constants::D2R;	  	  //linspace(start,stop,lamda_dstr). Goal is: symmetry				
			}
			//Uniformly dstr
			else if (argv[2]==string_uniform)   			
			{
				std::uniform_real_distribution <real> uniform_lamda(lamda_start_d, lamda_end_d); //All latitudes in this range are equally probable. Goal is: randomness
				lamda0  = uniform_lamda(generator) * Constants::D2R; 	
			}
			//Normaly dstr
			else if (argv[2]==string_normal)   			
			{
				real number;
				do
				{
					//Adaptive stdev according to the lamda domain range. Domain range (min,max) = (0,180).
					//That way for bigger domain ranges we have extended gausian, while for les wide ranges we have more "sharp".
					real stdev = (std::abs(lamda_end_d - lamda_start_d) / 180 ) * Constants::stdev_lamda;
					std::normal_distribution <real> normal_lamda(Constants::mean_lamda, stdev);
					number = normal_lamda(generator);

				}while(number<lamda_start_d || number>lamda_end_d); //Until valid lamda=(0,180).
				lamda0  = number * Constants::D2R; 	
			}
			else if (argv[2]==string_test)
			{
				lamda0 = Constants::lamda0;
			}
			else return EXIT_FAILURE;
		//--------------Latitude dstr--------------//


			const real Beq0 = Bmag_dipole(0);   		    //Beq isn't always Beq0?
			Blam0           = Bmag_dipole(lamda0);
			salpha0         = sin(aeq0)*sqrt(Blam0/Beq0);  //salpha = sin(aeq)*sqrt(Blam/Beq)
			if(aeq0*Constants::R2D>90) k=1;				   //Both k=1 and k=0 are valid for every particle!!(?). This basically defines if its upward or downward.
			else					   k=0;				   //We do this to distribute them, half upwards, half downwards.
			alpha0 = pow(-1,k)*asin(salpha0)+k*M_PI;       // sinx = a => x=(-1)^k * asin(a) + k*pi
			//Initialize and print particle state for this equatorial P.A and latitude.
			dstr[p].initialize(Constants::eta0,aeq0,alpha0,lamda0,Constants::Ekev0,Blam0,0,0,lamda_start_d,lamda_end_d);
			//std::cout<<"\nParticle"<<p<<" aeq0: "<< aeq0*Constants::R2D <<", lamda0: "<< lamda0*Constants::R2D;	
			p++; lamda_count++; //Next particle

		}	

		aeq_count++;
	}
	std::cout<<"\nDistribution done!";
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//
	//AEQ0 DISTRIBUTION CHECK
	const int sector_range = 10;
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
	
	std::vector<real> lamda_dstr(Constants::population), alpha_dstr(Constants::population), aeq_dstr(Constants::population), upar_dstr(Constants::population), uper_dstr(Constants::population), ppar_dstr(Constants::population), pper_dstr(Constants::population), eta_dstr(Constants::population), M_adiabatic_dstr(Constants::population), time_dstr(Constants::population), Ekin_dstr(Constants::population), zeta_dstr(Constants::population), trapped_dstr(Constants::population), escaped_dstr(Constants::population);

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
		trapped_dstr[p]    = dstr[p].trapped;
		escaped_dstr[p]    = dstr[p].escaped;
	}
	
	std::string file_name = "h5files/" + std::to_string(Constants::population) + "p_" + std::string(argv[1]) +"AEQ_" + std::string(argv[2]) + "LAMDA.h5";
	h5::File file(file_name, h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

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
	h5::DataSet data_trapped        = file.createDataSet("trapped", trapped_dstr);
	h5::DataSet data_escaped        = file.createDataSet("escaped", escaped_dstr);
	h5::DataSet aeq0bins            = file.createDataSet("aeq0_bins", aeq0_bins);

//----------------------------------------WRITE TO HDF5 FILE------------------------------------//

    return 0;
}
