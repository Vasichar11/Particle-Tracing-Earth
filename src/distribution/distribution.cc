#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <random> 
#include <iomanip>  
#include <omp.h>
#include <H5Cpp.h>
#include <algorithm>

#include "common.h"
#include "../headers/struct_Particles.h"
#include "../headers/functions.h"
#include "constants.h"


//Member function to find valid latitude range
void latitude_domain(real aeq0, real &latitude_start_d, real &latitude_end_d)
{
	const real Beq0 = Bmag_dipole(0);   	 //Beq isn't always Beq0?
	latitude_end_d     = 0;
	latitude_start_d   = 0;
	real salpha0;
	do
	{
		latitude_end_d += Lat_dstr::domain_step; 
		//Gradually increase and check if salpha0 is valid.
		//("Brute Force domain validation")
		real Blam0   = Bmag_dipole(latitude_end_d*Universal::D2R);
		salpha0      = sin(aeq0)*sqrt(Blam0/Beq0);  //salpha = sin(aeq)*sqrt(Blam/Beq)
							  
	}					//e.g aeq0_deg=89.9989367479725, latitude0_deg=0.000514656086141539 
	while( (salpha0<=1) && (salpha0>=-1) && (salpha0!=0) ) ; // &&latitude0<M_PI
	//Else Invalid
	/*-1<salpha || salpha>1 => domain error
	aeq0=0||180 => salpha0 = 0 => alpha0 = 0 => pper0 = 0 => PROBLEM IN Bell_param, a2(division with 0).
	sin(pi) is actually not zero in C++... */
	latitude_end_d   -= Lat_dstr::domain_step;     //This is the last (positive) valid value.
	latitude_start_d  = -latitude_end_d;	 //This would be the first (negative) valid value. 		
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

	std::cout<<"\n\nParticle population: " << Distribution::population;
	std::cout<<"\nDistributed as: |From    To|"<<"    (evenly, uniform, normal, test) \n";
	std::cout<<"\nEquatorial P.A: |"<<Aeq_dstr::start_deg<<" "<<Aeq_dstr::end_deg<<"|      "<< argv[1];
	std::cout<<"\nLatitude:       |"<<Lat_dstr::start_deg<<" "<<Lat_dstr::end_deg<<"|     "<< argv[2];
	std::cout<<"\nEta:            |"<<Eta_dstr::start_deg<<" "<<Eta_dstr::end_deg<<"|      "<< argv[3];
	std::cout<<"\nEkin:           |"<<Ekin_dstr::start<<" "<<Ekin_dstr::end<<"|    "<< argv[4]<< std::endl;



	Particles single; //Single particle struct.
	std::vector<Particles> dstr(Distribution::population, single);	//Vector of structs for particle distribution.
	real latitude0,aeq0,eta0,Ekin0;
	latitude0=aeq0=eta0=Ekin0=0;
	//const real Beq0 = Bmag_dipole(0);   	 //Beq isn't always Beq0?
	std::random_device seed;         //Random seed. 
	std::mt19937 generator(seed());  //PRNG initialized with seed.
	int p=0;
	int aeq_count = 0 ;
	int latitude_count = 0 ;
	int eta_count = 0;
	int Ekin_count = 0;


	for(int particle_count=0; particle_count<Distribution::population; particle_count++)
	{
		//------------------------------------------- DISTRIBUTE P.A ------------------------------------------------//
		// Evenly with goal: symmetry	
		if(argv[1]==string_evenly)   
		{
			// Step for linspace(start, end, aeq_groups) 
			const real aeq_step_d	  = (Aeq_dstr::end_deg - Aeq_dstr::start_deg)/(Aeq_dstr::steps-1); 	
			aeq0 = (Aeq_dstr::start_deg + aeq_count*aeq_step_d) * Universal::D2R;		
			aeq_count++;	
		}
		// Uniformly with goal: randomness
		else if (argv[1]==string_uniform)   			
		{
			//All values in this range are equally probable.
			std::uniform_real_distribution <real> uniform_aeq  (Aeq_dstr::start_deg, Aeq_dstr::end_deg);
			aeq0  = uniform_aeq(generator) * Universal::D2R; 	
		}
		// Normally with goal: reach normal dstr in logarithm scale
		else if (argv[1]==string_normal)   			
		{
			real number;
			do
			{
				std::normal_distribution <real> normal_aeq(Aeq_dstr::mean, Aeq_dstr::stdev);
				number = normal_aeq(generator);

			}while(number<=0 || number>=180); //Until valid aeq=(0,180).
			aeq0  = number * Universal::D2R; 	
		}
		// Constant with goal: testing
		else if (argv[1]==string_constant)
		{
			aeq0 = Aeq_dstr::value;
		}
		//----------------------------------------- DISTRIBUTE LATITUDE ----------------------------------------------//
		//-------------Latitude domain-------------//
		real latitude_end_d, latitude_start_d;
		latitude_domain(aeq0, latitude_start_d, latitude_end_d); //Finds latitude_start_d && latitude_end_d. 
		//std::cout<<"\nFor aeq0 = "<<aeq0*Universal::R2D<<" degrees\nThe latitude domain in degrees"<<"\n|From "<<" To|\n| "<<latitude_start_d << "  "<< " " << latitude_end_d <<"|\n";
		//-------------Latitude domain-------------//
		// Evenly with goal: symmetry	
		if(argv[2]==string_evenly)   
		{
			//Step for linspace(start,end,aeq_groups) 
			const real latitude_step_d	  = (latitude_end_d - latitude_start_d)/(Lat_dstr::steps-1);						
			latitude0 = (latitude_start_d + latitude_count*latitude_step_d) * Universal::D2R;	  					
		}
		// Uniformly with goal: randomness
		else if (argv[2]==string_uniform)   			
		{
			//All latitudes in this range are equally probable. 
			std::uniform_real_distribution <real> uniform_latitude(latitude_start_d, latitude_end_d); 
			latitude0  = uniform_latitude(generator) * Universal::D2R; 							 
		}
		// Normally with goal: reach normal dstr in logarithm scale.
		else if (argv[2]==string_normal)   			
		{
			real number;
			do
			{
				//Adaptive stdev according to the latitude domain range. Domain range (min,max) = (0,90+90).
				//That way for bigger domain ranges we have extended gausian, while for less wide ranges we have more "sharp".
				real stdev =((latitude_end_d-latitude_start_d)/180) * Lat_dstr::max_stdev;
				std::normal_distribution <real> normal_latitude(Lat_dstr::mean, stdev);
				number = normal_latitude(generator);

			}while(number<latitude_start_d || number>latitude_end_d); //Until valid latitude=(-90,90).
			latitude0  = number * Universal::D2R; 	
		}
		// Constant with goal: testing
		else if (argv[2]==string_constant)
		{
			latitude0 = Lat_dstr::value;
		}
		//-------------------------------------------- DISTRIBUTE ETA --------------------------------------------------//
		// Evenly with goal: symmetry	
		if(argv[3]==string_evenly)   
		{
			//Step for linspace(start,end,eta_steps) 
			const real eta_step_d	  = (Eta_dstr::end_deg - Eta_dstr::start_deg)/(Eta_dstr::steps-1); 	
			eta0 = (Eta_dstr::start_deg + eta_count*eta_step_d) * Universal::D2R;			
		}
		// Uniformly with goal: randomness
		else if (argv[3]==string_uniform)   			
		{
			//All values in this range are equally probable.
			std::uniform_real_distribution <real> uniform_eta  (Eta_dstr::start_deg, Eta_dstr::end_deg);
			eta0  = uniform_eta(generator) * Universal::D2R; 	
		}
		// Normally with goal: reach normal dstr 
		else if (argv[3]==string_normal)   			
		{
			real number;
			do
			{
				std::normal_distribution <real> normal_eta(Eta_dstr::mean, Eta_dstr::stdev);
				number = normal_eta(generator);

			}while(number<=0 || number>=180); //Until valid eta=(0,180).
			eta0  = number * Universal::D2R; 	
		}
		// Constant with goal: testing
		else if (argv[3]==string_constant)
		{
			eta0 = Eta_dstr::value;
		}
		//-------------------------------------------- DISTRIBUTE EKIN -------------------------------------------------//
		// Evenly with goal: symmetry	
		if(argv[4]==string_evenly)   
		{
			//Step for linspace(start,end,Ekin_steps) 
			const real Ekin_step = (Ekin_dstr::end - Ekin_dstr::start)/(Ekin_dstr::steps-1); 	
			Ekin0 = (Ekin_dstr::start + Ekin_count*Ekin_step);			
		}
		// Uniformly with goal: randomness
		else if (argv[4]==string_uniform)   			
		{
			//All values in this range are equally probable.
			std::uniform_real_distribution <real> uniform_Ekin  (Ekin_dstr::start, Ekin_dstr::end);
			Ekin0  = uniform_Ekin(generator); 	
		}
		// Normally with goal: reach normal dstr 
		else if (argv[4]==string_normal)   			
		{
			real number;
			do
			{
				std::normal_distribution <real> normal_Ekin(Ekin_dstr::mean, Ekin_dstr::stdev);
				number = normal_Ekin(generator);

			}while(number<=0 || number>=180); //Until valid Ekin=(0,180).
			Ekin0  = number; 	
		}
		// Constant with goal: testing
		else if (argv[4]==string_constant)
		{
			Ekin0 = Ekin_dstr::value;
		}

		//-------------------------------------------- FINALLY INITIALIZE -------------------------------------------------//

		dstr[p].initialize(eta0,aeq0,latitude0,Ekin0,0,0);
		// std::cout<<"\nParticle"<<p<<" aeq0: "<< dstr[p].aeq0*Universal::R2D <<", latitude0: "<< dstr[p].latitude0*Universal::R2D<<" eta0: "<< dstr[p].eta0*Universal::R2D <<", Ekin0: "<< dstr[p].Ekin0;	
		p++; //Next particle

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
	for(int p=0; p<Distribution::population; p++)
	{
		sec = floor((dstr[p].aeq0*Universal::R2D)/sector_range); //Which sector has this particle?
		aeq0_bins.at(sec) ++; 										  //This sector has this particle
	}
	std::cout<<"\nEquatorial P.A Initialization: ";
	for(int sector=0;sector<sectors;sector++)					      //Print population of P.A bins
	{				
		std::cout<<"\naeq0 range: "<<sector*sector_range<< " - " <<(sector+1)*sector_range<< " has " << aeq0_bins.at(sector) << " particles.";
	}

//---------------------------------------------------------------WRITE TO HDF5 FILE---------------------------------------------------------------------//
	// Assign from struct to 1d vectors. Trying std::execution::par

	// Create vectors
	std::vector<real> alpha0_dstr(Distribution::population);
	std::vector<real> latitude0_dstr(Distribution::population);
	std::vector<real> aeq0_dstr(Distribution::population);
	std::vector<real> ppar0_dstr(Distribution::population);
	std::vector<real> pper0_dstr(Distribution::population);
	std::vector<real> upar0_dstr(Distribution::population);
	std::vector<real> uper0_dstr(Distribution::population);
	std::vector<real> eta0_dstr(Distribution::population);
	std::vector<real> zeta0_dstr(Distribution::population);
	std::vector<real> Ekin0_dstr(Distribution::population);
	std::vector<real> M_adiabatic0_dstr(Distribution::population);
	std::vector<real> time0_dstr(Distribution::population);
	std::vector<real> trapped_dstr(Distribution::population);
	std::vector<real> escaped_dstr(Distribution::population);
	std::vector<real> negative_dstr(Distribution::population);
	std::vector<real> nan_dstr(Distribution::population);
	std::vector<real> high_dstr(Distribution::population);


	std::transform(dstr.begin(), dstr.end(), alpha0_dstr.begin(), [](const auto& d) { return d.alpha0; });
	std::transform(dstr.begin(), dstr.end(), latitude0_dstr.begin(), [](const auto& d) { return d.latitude0; });
	std::transform(dstr.begin(), dstr.end(), aeq0_dstr.begin(), [](const auto& d) { return d.aeq0; });
	std::transform(dstr.begin(), dstr.end(), ppar0_dstr.begin(), [](const auto& d) { return d.ppar0; });
	std::transform(dstr.begin(), dstr.end(), pper0_dstr.begin(), [](const auto& d) { return d.pper0; });
	std::transform(dstr.begin(), dstr.end(), upar0_dstr.begin(), [](const auto& d) { return d.upar0; });
	std::transform(dstr.begin(), dstr.end(), uper0_dstr.begin(), [](const auto& d) { return d.uper0; });
	std::transform(dstr.begin(), dstr.end(), eta0_dstr.begin(), [](const auto& d) { return d.eta0; });
	std::transform(dstr.begin(), dstr.end(), zeta0_dstr.begin(), [](const auto& d) { return d.zeta0; });
	std::transform(dstr.begin(), dstr.end(), Ekin0_dstr.begin(), [](const auto& d) { return d.Ekin0; });
	std::transform(dstr.begin(), dstr.end(), M_adiabatic0_dstr.begin(), [](const auto& d) { return d.M_adiabatic0; });
	std::transform(dstr.begin(), dstr.end(), time0_dstr.begin(), [](const auto& d) { return d.time0; });
	std::transform(dstr.begin(), dstr.end(), trapped_dstr.begin(), [](const auto& d) { return d.trapped; });
	std::transform(dstr.begin(), dstr.end(), escaped_dstr.begin(), [](const auto& d) { return d.escaped; });
	std::transform(dstr.begin(), dstr.end(), negative_dstr.begin(), [](const auto& d) { return d.negative; });
	std::transform(dstr.begin(), dstr.end(), nan_dstr.begin(), [](const auto& d) { return d.nan; });
	std::transform(dstr.begin(), dstr.end(), high_dstr.begin(), [](const auto& d) { return d.high; });

	
	// File name
	std::string file_name = "output/files/dstr" + std::to_string(Distribution::population) + "p_" + std::string(argv[1]) +"AEQ_" + std::string(argv[2]) + "LAMDA_" +std::string(argv[3])+"ETA_"+std::string(argv[4]) + "EKIN.h5";
	// Create instance of hdf5 file.
	//H5F_ACC_TRUNC flag -> new file or overwritei
	H5::H5File file(file_name, H5F_ACC_TRUNC);

	// Create dataspace and datatype
	hsize_t dim = Distribution::population;
	H5::DataSpace data_space(1,&dim);
	H5::DataType data_type(H5::PredType::NATIVE_DOUBLE);

	// Create the datasets
	H5::DataSet dataset_latitude0 = file.createDataSet("latitude0", data_type, data_space);
	H5::DataSet dataset_aeq0 = file.createDataSet("aeq0", data_type, data_space);
	H5::DataSet dataset_alpha0 = file.createDataSet("alpha0", data_type, data_space);
	H5::DataSet dataset_upar0 = file.createDataSet("upar0", data_type, data_space);
	H5::DataSet dataset_uper0 = file.createDataSet("uper0", data_type, data_space);
	H5::DataSet dataset_ppar0 = file.createDataSet("ppar0", data_type, data_space);
	H5::DataSet dataset_pper0 = file.createDataSet("pper0", data_type, data_space);
	H5::DataSet dataset_eta0 = file.createDataSet("eta0", data_type, data_space);
	H5::DataSet dataset_zeta0 = file.createDataSet("zeta0", data_type, data_space);
	H5::DataSet dataset_time0 = file.createDataSet("time0", data_type, data_space);
	H5::DataSet dataset_M_adiabatic0 = file.createDataSet("M_adiabatic0", data_type, data_space);
	H5::DataSet dataset_Ekin0 = file.createDataSet("Ekin0", data_type, data_space);
	H5::DataSet dataset_trapped0 = file.createDataSet("trapped0", data_type, data_space);
	H5::DataSet dataset_escaped0 = file.createDataSet("escaped0", data_type, data_space);
	H5::DataSet dataset_nan0 = file.createDataSet("nan0", data_type, data_space);
	H5::DataSet dataset_negative0 = file.createDataSet("negative0", data_type, data_space);
	H5::DataSet dataset_high0 = file.createDataSet("high0", data_type, data_space);
	H5::DataSet dataset_aeq0_bins = file.createDataSet("aeq0_bins", data_type, data_space);

	// Write
	dataset_latitude0.write(latitude0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_aeq0.write(aeq0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_alpha0.write(alpha0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_upar0.write(upar0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_uper0.write(uper0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_ppar0.write(ppar0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_pper0.write(pper0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_eta0.write(eta0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_zeta0.write(zeta0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_time0.write(time0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_M_adiabatic0.write(M_adiabatic0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_Ekin0.write(Ekin0_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_trapped0.write(trapped_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_escaped0.write(escaped_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_nan0.write(nan_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_negative0.write(negative_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	dataset_high0.write(high_dstr.data(), H5::PredType::NATIVE_DOUBLE);
	

	// Close datasets
	dataset_latitude0.close();
	dataset_aeq0.close();
	dataset_alpha0.close();
	dataset_upar0.close();
	dataset_uper0.close();
	dataset_ppar0.close();
	dataset_pper0.close();
	dataset_eta0.close();
	dataset_zeta0.close();
	dataset_time0.close();
	dataset_M_adiabatic0.close();
	dataset_Ekin0.close();
	dataset_trapped0.close();
	dataset_escaped0.close();
	dataset_nan0.close();
	dataset_negative0.close();
	dataset_high0.close();
	dataset_aeq0_bins.close();
	
	
	
//---------------------------------------------------------------WRITE TO HDF5 FILE: END---------------------------------------------------------------//


    return 0;
}
