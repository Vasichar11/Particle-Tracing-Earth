#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <random> 
#include <iomanip>  //For std::setprecision()
#include <omp.h>

#include "common.h"
#include "struct_Particles.h"   
#include "functions.h"
#include "constants.h"

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
		lamda_end_d += Constants::lamda_domain_step; 
		//Gradually increase and check if salpha0 is valid.
		//("Brute Force domain validation")
		real Blam0   = Bmag_dipole(lamda_end_d*Constants::D2R);
		salpha0      = sin(aeq0)*sqrt(Blam0/Beq0);  //salpha = sin(aeq)*sqrt(Blam/Beq)
							  
	}					//e.g aeq0_deg=89.9989367479725, lamda0_deg=0.000514656086141539 
	while( (salpha0<=1) && (salpha0>=-1) && (salpha0!=0) ) ; // &&lamda0<M_PI
	//Else Invalid
	/*-1<salpha || salpha>1 => domain error
	aeq0=0||180 => salpha0 = 0 => alpha0 = 0 => pper0 = 0 => PROBLEM IN Bell_param, a2(division with 0).
	sin(pi) is actually not zero in C++... */
	lamda_end_d   -= Constants::lamda_domain_step;     //This is the last (positive) valid value.
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
	std::string string_constant = "constant";
	
	//---ARGV ERROR---//
	if(argc!=5) throw std::invalid_argument("\nSet command line arguments from the list:(evenly, uniform, normal, test) to distribute the particles.\nargv[1] is for P.A, argv[2] is for latitude, arg[3] is for eta, argv[4] is for Ekin");
	for(int arg=1;arg<=4;arg++)
	{
		if( (argv[arg]!=string_evenly) && (argv[arg]!=string_uniform) && (argv[arg]!=string_normal) && (argv[arg]!=string_constant) ) 
			throw std::invalid_argument("\nArgument variables don't match any of the program's possible implementations.\nSet 4 command line arguments from the list:(evenly, uniform, normal, test) to distribute the particles.\nargv[1] for P.A, argv[2] for latitude, arg[3] for eta, argv[4] for Ekin");
	}

	std::cout<<"\n\nParticle population: " << Constants::population;
	std::cout<<"\nDistributed as: |From    To|"<<"    (evenly, uniform, normal, test) \n";
	std::cout<<"\nEquatorial P.A: |"<<Constants::aeq_start_d<<" "<<Constants::aeq_end_d<<"|      "<< argv[1];
	std::cout<<"\nLatitude:       |"<<Constants::lamda_start_d<<" "<<Constants::lamda_end_d<<"|     "<< argv[2];
	std::cout<<"\nEta:            |"<<Constants::eta_start_d<<" "<<Constants::eta_end_d<<"|      "<< argv[3];
	std::cout<<"\nEkin:           |"<<Constants::Ekin_start<<" "<<Constants::Ekin_end<<"|    "<< argv[4]<< std::endl;



	Particles single; //Single particle struct.
	std::vector<Particles> dstr(Constants::population, single);	//Vector of structs for particle distribution.
	real lamda0,aeq0,eta0,Ekin0;
	lamda0=aeq0=eta0=Ekin0=0;
	//const real Beq0 = Bmag_dipole(0);   	 //Beq isn't always Beq0?
	std::random_device seed;         //Random seed. 
	std::mt19937 generator(seed());  //PRNG initialized with seed.
	int p=0;
	


	for(int aeq_count=0; aeq_count<Constants::aeq_dstr; aeq_count++)
	{

		//------------------------------------------- DISTRIBUTE P.A ------------------------------------------------//
		//EVENLY. Goal is: symmetry	
		if(argv[1]==string_evenly)   
		{
			//Step for linspace(start,end,aeq_dstr) 
			const real aeq_step_d	  = (Constants::aeq_end_d - Constants::aeq_start_d)/(Constants::aeq_dstr-1); 	
			aeq0 = (Constants::aeq_start_d + aeq_count*aeq_step_d) * Constants::D2R;			
		}
		//UNIFORMLY. Goal is: randomness
		else if (argv[1]==string_uniform)   			
		{
			//All values in this range are equally probable.
			std::uniform_real_distribution <real> uniform_aeq  (Constants::aeq_start_d, Constants::aeq_end_d);
			aeq0  = uniform_aeq(generator) * Constants::D2R; 	
		}
		//NORMALLY. Goal is: reach normal dstr in logarithm scale.
		else if (argv[1]==string_normal)   			
		{
			real number;
			do
			{
				std::normal_distribution <real> normal_aeq(Constants::mean_aeq, Constants::stdev_aeq);
				number = normal_aeq(generator);

			}while(number<=0 || number>=180); //Until valid aeq=(0,180).
			aeq0  = number * Constants::D2R; 	
		}
		//CONSTANT. Goal is: testing
		else if (argv[1]==string_constant)
		{
			aeq0 = Constants::aeq0;
		}
		//------------------------------------------- DISTRIBUTE P.A ------------------------------------------------//


		//----------------------------------------- DISTRIBUTE LATITUDE ----------------------------------------------//
		
		//-------------Latitude domain-------------//
		real lamda_end_d, lamda_start_d;
		lamda_domain(aeq0, lamda_start_d, lamda_end_d); //Finds lamda_start_d && lamda_end_d. 
		//std::cout<<"\nFor aeq0 = "<<aeq0*Constants::R2D<<" degrees\nThe latitude domain in degrees"<<"\n|From "<<" To|\n| "<<lamda_start_d << "  "<< " " << lamda_end_d <<"|\n";
		//-------------Latitude domain-------------//


		for(int lamda_count=0; lamda_count<Constants::lamda_dstr; lamda_count++)
		{
			//EVENLY. Goal is: symmetry	
			if(argv[2]==string_evenly)   
			{
				//Step for linspace(start,end,aeq_dstr) 
				const real lamda_step_d	  = (lamda_end_d - lamda_start_d)/(Constants::lamda_dstr-1);						
				lamda0 = (lamda_start_d + lamda_count*lamda_step_d) * Constants::D2R;	  					
			}
			//UNIFORMLY. Goal is: randomness
			else if (argv[2]==string_uniform)   			
			{
				//All latitudes in this range are equally probable. 
				std::uniform_real_distribution <real> uniform_lamda(lamda_start_d, lamda_end_d); 
				lamda0  = uniform_lamda(generator) * Constants::D2R; 							 
			}
			//NORMALLY. Goal is: reach normal dstr in logarithm scale.
			else if (argv[2]==string_normal)   			
			{
				real number;
				do
				{
					//Adaptive stdev according to the lamda domain range. Domain range (min,max) = (0,90+90).
					//That way for bigger domain ranges we have extended gausian, while for less wide ranges we have more "sharp".
					real stdev =((lamda_end_d-lamda_start_d)/180) * Constants::max_stdev_lamda;
					std::normal_distribution <real> normal_lamda(Constants::mean_lamda, stdev);
					number = normal_lamda(generator);

				}while(number<lamda_start_d || number>lamda_end_d); //Until valid lamda=(-90,90).
				lamda0  = number * Constants::D2R; 	
			}
			//CONSTANT. Goal is: testing
			else if (argv[2]==string_constant)
			{
				lamda0 = Constants::lamda0;
			}
		//----------------------------------------- DISTRIBUTE LATITUDE ----------------------------------------------//

			
		//-------------------------------------------- DISTRIBUTE ETA --------------------------------------------------//
			for(int eta_count=0; eta_count<Constants::eta_dstr; eta_count++)
			{
				//EVENLY. Goal is: symmetry	
				if(argv[3]==string_evenly)   
				{
					//Step for linspace(start,end,eta_dstr) 
					const real eta_step_d	  = (Constants::eta_end_d - Constants::eta_start_d)/(Constants::eta_dstr-1); 	
					eta0 = (Constants::eta_start_d + eta_count*eta_step_d) * Constants::D2R;			
				}
				//UNIFORMLY. Goal is: randomness
				else if (argv[3]==string_uniform)   			
				{
					//All values in this range are equally probable.
					std::uniform_real_distribution <real> uniform_eta  (Constants::eta_start_d, Constants::eta_end_d);
					eta0  = uniform_eta(generator) * Constants::D2R; 	
				}
				//NORMALLY. Goal is: reach normal dstr 
				else if (argv[3]==string_normal)   			
				{
					real number;
					do
					{
						std::normal_distribution <real> normal_eta(Constants::mean_eta, Constants::stdev_eta);
						number = normal_eta(generator);

					}while(number<=0 || number>=180); //Until valid eta=(0,180).
					eta0  = number * Constants::D2R; 	
				}
				//CONSTANT. Goal is: testing
				else if (argv[3]==string_constant)
				{
					eta0 = Constants::eta0;
				}
		//-------------------------------------------- DISTRIBUTE ETA --------------------------------------------------//
				

		//-------------------------------------------- DISTRIBUTE EKIN -------------------------------------------------//
				for(int Ekin_count=0; Ekin_count<Constants::Ekin_dstr; Ekin_count++)
				{
					//EVENLY. Goal is: symmetry	
					if(argv[4]==string_evenly)   
					{
						//Step for linspace(start,end,Ekin_dstr) 
						const real Ekin_step = (Constants::Ekin_end - Constants::Ekin_start)/(Constants::Ekin_dstr-1); 	
						Ekin0 = (Constants::Ekin_start + Ekin_count*Ekin_step);			
					}
					//UNIFORMLY. Goal is: randomness
					else if (argv[4]==string_uniform)   			
					{
						//All values in this range are equally probable.
						std::uniform_real_distribution <real> uniform_Ekin  (Constants::Ekin_start, Constants::Ekin_end);
						Ekin0  = uniform_Ekin(generator); 	
					}
					//NORMALLY. Goal is: reach normal dstr 
					else if (argv[4]==string_normal)   			
					{
						real number;
						do
						{
							std::normal_distribution <real> normal_Ekin(Constants::mean_Ekin, Constants::stdev_Ekin);
							number = normal_Ekin(generator);

						}while(number<=0 || number>=180); //Until valid Ekin=(0,180).
						Ekin0  = number; 	
					}
					//CONSTANT. Goal is: testing
					else if (argv[4]==string_constant)
					{
						Ekin0 = Constants::Ekin0;
					}
		//-------------------------------------------- DISTRIBUTE EKIN -------------------------------------------------//
			

					//Initialize particle.
					dstr[p].initialize(eta0,aeq0,lamda0,Ekin0,0,0);
					std::cout<<"\nParticle"<<p<<" aeq0: "<< dstr[p].aeq0*Constants::R2D <<", lamda0: "<< dstr[p].lamda0*Constants::R2D<<" eta0: "<< dstr[p].eta0*Constants::R2D <<", Ekin0: "<< dstr[p].Ekin0;	
					p++; //Next particle

				}
			}
		}	
	}
	std::cout<<"\nDistribution done!";
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//
	//AEQ0 DISTRIBUTION CHECK
	const int sector_range = 1;
	const int view = 180;
	const int sectors = view/sector_range;
	std::array<int, sectors> aeq0_bins;
	std::fill(std::begin(aeq0_bins), std::end(aeq0_bins), 0); 		  //Initialize array elements with 0
	int sec;
	for(int p=0; p<Constants::population; p++)
	{
		sec = floor((dstr[p].aeq0*Constants::R2D)/sector_range); //Which sector has this particle?
		aeq0_bins.at(sec) ++; 										  //This sector has this particle
	}
	std::cout<<"\nEquatorial P.A Initialization: ";
	for(int sector=0;sector<sectors;sector++)					      //Print population of P.A bins
	{				
		std::cout<<"\naeq0 range: "<<sector*sector_range<< " - " <<(sector+1)*sector_range<< " has " << aeq0_bins.at(sector) << " particles.";
	}

//----------------------------------------WRITE TO HDF5 FILE------------------------------------//
	
	std::vector<real> lamda_dstr(Constants::population), alpha_dstr(Constants::population), aeq_dstr(Constants::population), upar_dstr(Constants::population), uper_dstr(Constants::population), ppar_dstr(Constants::population), pper_dstr(Constants::population), eta_dstr(Constants::population), M_adiabatic_dstr(Constants::population), time_dstr(Constants::population), Ekin_dstr(Constants::population), zeta_dstr(Constants::population), trapped_dstr(Constants::population), escaped_dstr(Constants::population), nan_dstr(Constants::population), negative_dstr(Constants::population), high_dstr(Constants::population);

	//Assign from struct to 1d vectors.
	for(int p=0; p<Constants::population; p++)
	{
		alpha_dstr[p]      = dstr[p].alpha0;
		lamda_dstr[p]      = dstr[p].lamda0;
		aeq_dstr[p]        = dstr[p].aeq0;
		ppar_dstr[p]       = dstr[p].ppar0;
		pper_dstr[p]       = dstr[p].pper0;
		upar_dstr[p]       = dstr[p].upar0;
		uper_dstr[p]       = dstr[p].uper0;
		eta_dstr[p]        = dstr[p].eta0;
		zeta_dstr[p]       = dstr[p].zeta0;
		Ekin_dstr[p]       = dstr[p].Ekin0;
		M_adiabatic_dstr[p]= dstr[p].M_adiabatic0;
		time_dstr[p]       = dstr[p].time0;
		trapped_dstr[p]    = dstr[p].trapped;
		escaped_dstr[p]    = dstr[p].escaped;
		negative_dstr[p]   = dstr[p].negative;
		nan_dstr[p]    	   = dstr[p].nan;
		high_dstr[p]	   = dstr[p].high;
	}
	
	std::string file_name = "output/files/" + std::to_string(Constants::population) + "p_" + std::string(argv[1]) +"AEQ_" + std::string(argv[2]) + "LAMDA_" +std::string(argv[3])+"ETA_"+std::string(argv[4]) + "EKIN.h5";
	h5::File file(file_name, h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

	h5::DataSet data_lat            = file.createDataSet("lamda0", lamda_dstr);
	h5::DataSet data_aeq            = file.createDataSet("aeq0", aeq_dstr);
	h5::DataSet data_alpha          = file.createDataSet("alpha0", alpha_dstr);
	h5::DataSet data_upar           = file.createDataSet("upar0", upar_dstr);
	h5::DataSet data_uper           = file.createDataSet("uper0", uper_dstr);
	h5::DataSet data_ppar           = file.createDataSet("ppar0", ppar_dstr);
	h5::DataSet data_pper           = file.createDataSet("pper0", pper_dstr);
	h5::DataSet data_eta            = file.createDataSet("eta0",  eta_dstr);
	h5::DataSet data_zeta           = file.createDataSet("zeta0", zeta_dstr);
	h5::DataSet data_time           = file.createDataSet("time0", time_dstr);
	h5::DataSet data_M_adiabatic    = file.createDataSet("M_adiabatic0", M_adiabatic_dstr);
	h5::DataSet data_Ekin           = file.createDataSet("Ekin0", Ekin_dstr);
	h5::DataSet data_trapped        = file.createDataSet("trapped0", trapped_dstr);
	h5::DataSet data_escaped        = file.createDataSet("escaped0", escaped_dstr);
	h5::DataSet data_nan       	    = file.createDataSet("nan0", nan_dstr);
	h5::DataSet data_negative       = file.createDataSet("negative0", negative_dstr);
	h5::DataSet data_high      	    = file.createDataSet("high0", high_dstr);
	h5::DataSet aeq0bins            = file.createDataSet("aeq0_bins", aeq0_bins);

//----------------------------------------WRITE TO HDF5 FILE------------------------------------//

    return 0;
}
