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

#include "headers/functions.h"

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
		aeq0_mr = M_PI - aeq0;		 													  		  //Mirror aeq  														
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
		{																				   //Loop until <lamda_dstr> valid states for this aeq0.
			std::uniform_real_distribution <real> distribution(-interval/2, interval/2);   //Uniform distribution with varying interval.
			number  = distribution(generator);											   //When aeq~90, interval is small. Moving away from 90, interval increases.
			lamda0 	= number * Constants::D2R;											    
			lamda0_mr 	= - lamda0;											 			   //Mirror lamda
			std::cout<<"\n"<<lamda0;

			//Find P.A at lamda0.
			Blam0 	   = Bmag_dipole(lamda0);
			Blam0_mr   = Bmag_dipole(lamda0_mr);
			salpha0    = sin(aeq0)*sqrt(Blam0/Beq0); 										//(2.20) Bortnik thesis
			salpha0_mr    = sin(aeq0_mr)*sqrt(Blam0_mr/Beq0);  
			if(   !((salpha0>1) || (salpha0<-1) || (salpha0==0) || (salpha0_mr>1) || (salpha0_mr<-1) || (salpha0_mr==0) )   )  
			{																			 	//NOR of these should be true. Otherwise domain error.
				//Projecting aeq from alpha
				k       = ((aeq0*   Constants::R2D>90)    ? 1 : 0);     					//kEN...(here k=0 or 1 ?)
				k_mr    = ((aeq0_mr*Constants::R2D>90)    ? 1 : 0); 		
				alpha0  = pow(-1,k)*asin(salpha0)+k*M_PI;			 						// sinx = a => x=(-1)^k * asin(a) + k*pi
				alpha0_mr  = pow(-1,k_mr)*asin(salpha0_mr)+k_mr*M_PI;			 	
				dstr[p].initialize(Constants::eta0,aeq0,alpha0,lamda0,Constants::Ekev0,Blam0,0,0,0);
				dstr[p+1].initialize(Constants::eta0,aeq0_mr,alpha0_mr,lamda0_mr,Constants::Ekev0,Blam0_mr,0,0,0);
				//Print initial state of particles.
				//std::cout<<"\nParticle"<<p<<" aeq0: "<< aeq0*Constants::R2D <<", lamda0: "<< lamda0*Constants::R2D <<" gives alpha0: "<<alpha0*Constants::R2D<<std::endl;	
				//std::cout<<"\nParticle"<<p+1<<" aeq0: "<< aeq0_mr*Constants::R2D <<", lamda0: "<< lamda0_mr*Constants::R2D <<" gives alpha0: "<<alpha0_mr*Constants::R2D<<std::endl;
				p+=2;
				lat_count+=2;
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
	h5::File file("h5files/distribution_5000.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

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

#include "headers/read_csv.h"
#include "headers/interpolate.h"
#include "headers/common.h"
#include "headers/constants.h"
#include <vector>
#include <chrono>
#include <iomanip>  

#include <highfive/H5File.hpp>                                      
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;

int main()

{


    auto start1 = std::chrono::high_resolution_clock::now();

//-------------------------------------------------- READ CSV --------------------------------------------------------// 
    //Using function "rcsv" to read Ray data && store values in a single(!) vector.

    std::vector<std::pair<std::string, std::vector<real>>> ray_tracing = rcsv("L2_freq2500_psi-89_lat_0_damping.csv");

    int column_size = ray_tracing.at(0).second.size() ;       //Size of columns

    //Declaring seperate vectors to store Ray data individually.
    std::vector <real> timef, posx, posy, posz, vprelx, vprely, vprelz, vgrelx, vgrely, vgrelz, nx, ny, nz;
    std::vector <real> Bx, By, Bz, w, Ns1, Ns2, Ns3, Ns4, psi, theta_res, Y, Lf, alt, lat, lon, damp;
    std::vector <real> S_stix, D_stix, P_stix, R_stix, L_stix, wmega_e ;
    std::vector <real> Nspec,qs1,qs2,qs3,qs4,ms1,ms2,ms3,ms4,nus1,nus2,nus3,nus4;


    //Allocation of values inside seperate vectors 
    for(int i=0; i<column_size; i++)                
    {   
        timef.push_back(ray_tracing.at(0).second.at(i));
        posx.push_back(ray_tracing.at(1).second.at(i));       //Position of wave
        posy.push_back(ray_tracing.at(2).second.at(i));       //
        posz.push_back(ray_tracing.at(3).second.at(i));       //
        vprelx.push_back(ray_tracing.at(4).second.at(i));     //Phase velocity of wave
        vprely.push_back(ray_tracing.at(5).second.at(i));     //
        vprelz.push_back(ray_tracing.at(6).second.at(i));     //
        vgrelx.push_back(ray_tracing.at(7).second.at(i));     //Group velocity of wave
        vgrely.push_back(ray_tracing.at(8).second.at(i));     //
        vgrelz.push_back(ray_tracing.at(9).second.at(i));     //
        nx.push_back(ray_tracing.at(10).second.at(i));        //Refractive index of propagation medium
        ny.push_back(ray_tracing.at(11).second.at(i));        //
        nz.push_back(ray_tracing.at(12).second.at(i));        //
        Bx.push_back(ray_tracing.at(13).second.at(i));        //Magnetic field components of wave
        By.push_back(ray_tracing.at(14).second.at(i));        //
        Bz.push_back(ray_tracing.at(15).second.at(i));        //                                                                        
        w.push_back(ray_tracing.at(16).second.at(i));         //(constant) Angular frequency of wave
        Nspec.push_back(ray_tracing.at(17).second.at(i));     //(constant) Number of species, ions are single charged.
        qs1.push_back(ray_tracing.at(18).second.at(i));       //(constant) Charge of electron particle
        qs2.push_back(ray_tracing.at(19).second.at(i));       //(constant) Charge of positive Ions    
        qs3.push_back(ray_tracing.at(20).second.at(i));       //(constant)
        qs4.push_back(ray_tracing.at(21).second.at(i));       //(constant)
        ms1.push_back(ray_tracing.at(22).second.at(i));       //(constant) Mass of electron particle 
        ms2.push_back(ray_tracing.at(23).second.at(i));       //(constant) Mass of hydrogen particle 
        ms3.push_back(ray_tracing.at(24).second.at(i));       //(constant) Mass of helium particle    
        ms4.push_back(ray_tracing.at(25).second.at(i));       //(constant) Mass of oxygen particle    
        Ns1.push_back(ray_tracing.at(26).second.at(i));       //Density of electrons in radiation belt
        Ns2.push_back(ray_tracing.at(27).second.at(i));       //Density of hydrogen
        Ns3.push_back(ray_tracing.at(28).second.at(i));       //
        Ns4.push_back(ray_tracing.at(29).second.at(i));       //
        nus1.push_back(ray_tracing.at(30).second.at(i));      //(zero) Collision frequency of electron(?) particles
        nus2.push_back(ray_tracing.at(31).second.at(i));      //(zero) 
        nus3.push_back(ray_tracing.at(32).second.at(i));      //(zero)
        nus4.push_back(ray_tracing.at(33).second.at(i));      //(zero)
        psi.push_back(ray_tracing.at(34).second.at(i));       //Wave normal vector
        theta_res.push_back(ray_tracing.at(35).second.at(i)); //Resonance angle
        Y.push_back(ray_tracing.at(36).second.at(i));         //psi - theta_res
        Lf.push_back(ray_tracing.at(37).second.at(i));        //L_shell
        alt.push_back(ray_tracing.at(38).second.at(i));       //Altitude
        lat.push_back(ray_tracing.at(39).second.at(i));       //Latitude
        lon.push_back(ray_tracing.at(40).second.at(i));       //Longitude
        damp.push_back(ray_tracing.at(41).second.at(i));      //damping factor?
        S_stix.push_back(ray_tracing.at(42).second.at(i));    //Stix parameters
        D_stix.push_back(ray_tracing.at(43).second.at(i));    //
        P_stix.push_back(ray_tracing.at(44).second.at(i));    //
        R_stix.push_back(ray_tracing.at(45).second.at(i));    //
        L_stix.push_back(ray_tracing.at(46).second.at(i));    //
        wmega_e.push_back(ray_tracing.at(47).second.at(i));   //Gyrofrequency of electrons
    }                                                         //-----------Interpolation in constant vectors?---------//
//-------------------------------------------------- READ CSV: DONE --------------------------------------------------//
    auto stop1= std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    real time1 = duration1.count()*pow(10,-6);
    std::cout << "\nReading CSV. Time elapsed: " << time1 <<" seconds" ;

    auto start2 = std::chrono::high_resolution_clock::now();

//------------------------------------------------- INTERPOLATE RAY --------------------------------------------------//
    //std::cout << timef.at(0) << " " << timef.back();                  //Csv's first and last time value.
        
    int t_size = timef.size();      
    real last_timestep = timef.at(t_size - 2);                          //Second to last element.
    last_timestep = (int)(last_timestep*100 + 0.5) ;                    //Integer (value*100 + 0.5)
    last_timestep = last_timestep/100 ;                                 //---> Round last_timestep, 2 decimals.
    //std::cout<< "\nlast timestep " << last_timestep ;                 //Simulation's last timestep.
    
    std::vector <real> time_new {0} ;                                   //Initialize to use .back() in first iteration.
    int j=1;
    while(time_new.back()+(Constants::h) < last_timestep - Constants::h) //Arrange, 0 to last_timestep with stepsize h
    {                                //Why last_timestep - h ??
        time_new.push_back(0+j*Constants::h); 
        j++ ;
    }  

    std::vector <real> posx_int      =   interpolate(time_new, timef, posx);
    std::vector <real> posy_int      =   interpolate(time_new, timef, posy);
    std::vector <real> posz_int      =   interpolate(time_new, timef, posz);
    std::vector <real> vprelx_int    =   interpolate(time_new, timef, vprelx);
    std::vector <real> vprely_int    =   interpolate(time_new, timef, vprely);
    std::vector <real> vprelz_int    =   interpolate(time_new, timef, vprelz);
    std::vector <real> vgrelx_int    =   interpolate(time_new, timef, vgrelx);
    std::vector <real> vgrely_int    =   interpolate(time_new, timef, vgrely);
    std::vector <real> vgrelz_int    =   interpolate(time_new, timef, vgrelz);
    std::vector <real> nx_int        =   interpolate(time_new, timef, nx);
    std::vector <real> ny_int        =   interpolate(time_new, timef, ny);
    std::vector <real> nz_int        =   interpolate(time_new, timef, nz);
    std::vector <real> Bx_int        =   interpolate(time_new, timef, Bx);
    std::vector <real> By_int        =   interpolate(time_new, timef, By);
    std::vector <real> Bz_int        =   interpolate(time_new, timef, Bz);
    std::vector <real> w_int         =   interpolate(time_new, timef, w);
    std::vector <real> qs1_int       =   interpolate(time_new, timef, qs1);
    std::vector <real> qs2_int       =   interpolate(time_new, timef, qs2);
    std::vector <real> qs3_int       =   interpolate(time_new, timef, qs3);
    std::vector <real> qs4_int       =   interpolate(time_new, timef, qs4);
    std::vector <real> ms1_int       =   interpolate(time_new, timef, ms1);
    std::vector <real> ms2_int       =   interpolate(time_new, timef, ms2);
    std::vector <real> ms3_int       =   interpolate(time_new, timef, ms3);
    std::vector <real> ms4_int       =   interpolate(time_new, timef, ms4);
    std::vector <real> Ns1_int       =   interpolate(time_new, timef, Ns1);
    std::vector <real> Ns2_int       =   interpolate(time_new, timef, Ns2);
    std::vector <real> Ns3_int       =   interpolate(time_new, timef, Ns3);
    std::vector <real> Ns4_int       =   interpolate(time_new, timef, Ns4);
    std::vector <real> psi_int       =   interpolate(time_new, timef, psi);
    std::vector <real> theta_res_int =   interpolate(time_new, timef, theta_res);
    std::vector <real> Y_int         =   interpolate(time_new, timef, Y);
    std::vector <real> Lf_int        =   interpolate(time_new, timef, Lf);
    std::vector <real> alt_int       =   interpolate(time_new, timef, alt);
    std::vector <real> lat_int       =   interpolate(time_new, timef, lat);
    std::vector <real> lon_int       =   interpolate(time_new, timef, lon);
    std::vector <real> damp_int      =   interpolate(time_new, timef, damp);
    std::vector <real> S_stix_int    =   interpolate(time_new, timef, S_stix);
    std::vector <real> D_stix_int    =   interpolate(time_new, timef, D_stix);
    std::vector <real> P_stix_int    =   interpolate(time_new, timef, P_stix);
    std::vector <real> R_stix_int    =   interpolate(time_new, timef, R_stix);
    std::vector <real> L_stix_int    =   interpolate(time_new, timef, L_stix);
    std::vector <real> wmega_e_int   =   interpolate(time_new, timef, wmega_e);


//---------------------------------------------- INTERPOLATE RAY: DONE -----------------------------------------------//
    auto stop2= std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
    real time2 = duration2.count()*pow(10,-6);
    std::cout << "\nInterpolating Ray. Time elapsed: " << time2 <<" seconds" ;

    auto start3 = std::chrono::high_resolution_clock::now();

//------------------------------------------- OTHER VECTOR CALCULATIONS --------------------------------------------- //

  

    //Declare with known size
    std::vector <real> mu_ray, mu_sq, spsi, cpsi, kx_ray, kz_ray, kappa_ray, X_stix;
    std::vector <real> rho1, rho2, Byw_sq, fac1, Byw, Bxw, Bzw, Exw, Eyw, Ezw, Bw_ray, w1, w2, R1, R2 ;

for(size_t i=0; i<time_new.size(); i++)
{
    //Calculate more parameters #In[6]:
    mu_ray.push_back( sqrt(nx_int.at(i)*nx_int.at(i)+ny_int.at(i)*ny_int.at(i)+nz_int.at(i)*nz_int.at(i)) );
    mu_sq.push_back(  mu_ray.at(i)*mu_ray.at(i) );
    spsi.push_back( sin(psi_int.at(i)*Constants::D2R) );
    cpsi.push_back( cos(psi_int.at(i)*Constants::D2R) );
    kx_ray.push_back( (nx_int.at(i)*w_int.at(i)) / Constants::c );
    kz_ray.push_back( (nz_int.at(i)*w_int.at(i)) / Constants::c );
    kappa_ray.push_back( (mu_ray.at(i)*w_int.at(i)) / Constants::c );
    X_stix.push_back( P_stix_int.at(i)/(P_stix_int.at(i)-mu_sq.at(i)*spsi.at(i)*spsi.at(i)) );
    rho1.push_back( ((mu_sq.at(i)-S_stix_int.at(i))*mu_sq.at(i)*spsi.at(i)*cpsi.at(i)) / (D_stix_int.at(i)*(mu_sq.at(i)*spsi.at(i)*spsi.at(i)-P_stix_int.at(i))) );
    rho2.push_back( (mu_sq.at(i)-S_stix_int.at(i)) / D_stix_int.at(i) );

    //Define whistler waves #In[7]:R1
    Byw_sq.push_back( ((2.0*Constants::mu_0/Constants::c)*((Constants::pwr*damp_int.at(i))*X_stix.at(i)*X_stix.at(i)*rho2.at(i)*rho2.at(i)*std::abs(cos(psi_int.at(i)*Constants::D2R)))
                    /sqrt(pow((tan(psi_int.at(i)*Constants::D2R)-rho1.at(i)*rho2.at(i)*X_stix.at(i)),2) + pow((1+rho2.at(i)*rho2.at(i)*X_stix.at(i)),2))) );
    fac1.push_back( (P_stix_int.at(i)-mu_sq.at(i)*pow(sin(psi_int.at(i)*Constants::D2R),2)) );
    Byw.push_back( sqrt(Byw_sq.at(i)) );
    Bxw.push_back( std::abs((-(D_stix_int.at(i)*fac1.at(i))/(P_stix_int.at(i)*(S_stix_int.at(i)-mu_sq.at(i))))*Byw.at(i)) );
    Bzw.push_back( std::abs(((D_stix_int.at(i)*sin(psi_int.at(i)*Constants::D2R)*fac1.at(i))/(P_stix_int.at(i)*cos(psi_int.at(i)*Constants::D2R)*(S_stix_int.at(i)-mu_sq.at(i))))*Byw.at(i)) );
    Exw.push_back( std::abs(((Constants::c*fac1.at(i))/(mu_ray.at(i)*P_stix_int.at(i)*cos(psi_int.at(i)*Constants::D2R))*Byw.at(i))) );
    Eyw.push_back( std::abs(((D_stix_int.at(i)*Constants::c*fac1.at(i))/(mu_ray.at(i)*P_stix_int.at(i)*cos(psi_int.at(i)*Constants::D2R)*(pow(mu_ray.at(i),2)-S_stix_int.at(i))))*Byw.at(i)) );
    Ezw.push_back( std::abs((-(Constants::c*mu_ray.at(i)*sin(psi_int.at(i)))/P_stix_int.at(i))*Byw.at(i)) );

    Bw_ray.push_back( sqrt(Bxw.at(i)*Bxw.at(i) + Byw.at(i)*Bzw.at(i) + Byw.at(i)*Bzw.at(i)) );
    //From Bell parameters
    w1.push_back( (Constants::q_e/(2*Constants::m_e))*(Bxw.at(i)+Byw.at(i)) );
    w2.push_back( (Constants::q_e/(2*Constants::m_e))*(Bxw.at(i)-Byw.at(i)) );
    R1.push_back( (Exw.at(i)+Eyw.at(i))/(Bxw.at(i)+Byw.at(i)) );   
    R2.push_back( (Exw.at(i)-Eyw.at(i))/(Bxw.at(i)-Byw.at(i)) );
}
//--------------------------------------- OTHER VECTOR CALCULATIONS: DONE ------------------------------------------- //
    auto stop3 = std::chrono::high_resolution_clock::now();
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
    real time3 = duration3.count()*pow(10,-6);
    std::cout << "\nCalculating other vectors. Time elapsed: " << time3 <<" seconds" ;

    auto start4 = std::chrono::high_resolution_clock::now();


//------------------------------------------------- WRITE HDF5 FILE ------------------------------------------------- //

    h5::File file_out("h5files/interpolated_ray.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

    //Vectors from interpolation
    //h5::DataSet dataset_posx_int =      file_out.createDataSet("posx_int", posx_int);
    //h5::DataSet dataset_posy_int =      file_out.createDataSet("posy_int", posy_int);
    //h5::DataSet dataset_posz_int =      file_out.createDataSet("posz_int",posz_int);
    //h5::DataSet dataset_vprelx_int =    file_out.createDataSet("vprelx_int", vprelx_int);
    //h5::DataSet dataset_vprely_int =    file_out.createDataSet("vprely_int", vprely_int);
    //h5::DataSet dataset_vprelz_int =    file_out.createDataSet("vprelz_int", vprelz_int);
    //h5::DataSet dataset_vgrelx_int =    file_out.createDataSet("vgrelx_int", vgrelx_int);
    //h5::DataSet dataset_vgrely_int =    file_out.createDataSet("vgrely_int", vgrely_int);
    //h5::DataSet dataset_vgrelz_int =    file_out.createDataSet("vgrelz_int", vgrelz_int);
    //h5::DataSet dataset_nx_int =        file_out.createDataSet("nx_int", nx_int);
    //h5::DataSet dataset_ny_int =        file_out.createDataSet("ny_int", ny_int);
    //h5::DataSet dataset_nz_int =        file_out.createDataSet("nz_int", nz_int);
    //h5::DataSet dataset_Bx_int =        file_out.createDataSet("Bx_int", Bx_int);
    //h5::DataSet dataset_By_int =        file_out.createDataSet("By_int", By_int);
    //h5::DataSet dataset_Bz_int =        file_out.createDataSet("Bz_int", Bz_int);
    //h5::DataSet dataset_w_int =         file_out.createDataSet("w_int", w_int);
    //h5::DataSet dataset_qs1_int =       file_out.createDataSet("qs1_int", qs1_int);
    //h5::DataSet dataset_qs2_int =       file_out.createDataSet("qs2_int", qs2_int);
    //h5::DataSet dataset_qs3_int =       file_out.createDataSet("qs3_int", qs3_int);
    //h5::DataSet dataset_qs4_int =       file_out.createDataSet("qs4_int", qs4_int);
    //h5::DataSet dataset_ms1_int =       file_out.createDataSet("ms1_int", ms1_int);
    //h5::DataSet dataset_ms2_int =       file_out.createDataSet("ms2_int", ms2_int);
    //h5::DataSet dataset_ms3_int =       file_out.createDataSet("ms3_int", ms3_int);
    //h5::DataSet dataset_ms4_int =       file_out.createDataSet("ms4_int", ms4_int);
    //h5::DataSet dataset_Ns1_int =       file_out.createDataSet("Ns1_int", Ns1_int);
    //h5::DataSet dataset_Ns2_int =       file_out.createDataSet("Ns2_int", Ns2_int);
    //h5::DataSet dataset_Ns3_int =       file_out.createDataSet("Ns3_int", Ns3_int);
    //h5::DataSet dataset_Ns4_int =       file_out.createDataSet("Ns4_int", Ns4_int);
    //h5::DataSet dataset_psi_int =       file_out.createDataSet("psi_int", psi_int);
    //h5::DataSet dataset_theta_res_int = file_out.createDataSet("theta_res_int", theta_res_int);
    //h5::DataSet dataset_Y_int =         file_out.createDataSet("Y_int", Y_int);
    //h5::DataSet dataset_Lf_int =        file_out.createDataSet("Lf_int", Lf_int);
    //h5::DataSet dataset_alt_int =       file_out.createDataSet("alt_int", alt_int);
    h5::DataSet dataset_lat_int =       file_out.createDataSet("lat_int", lat_int);
    //h5::DataSet dataset_lon_int =       file_out.createDataSet("lon_int", lon_int);
    //h5::DataSet dataset_damp_int =      file_out.createDataSet("damp_int", damp_int);
    //h5::DataSet dataset_S_stix_int =    file_out.createDataSet("S_stix_int", S_stix_int);
    //h5::DataSet dataset_D_stix_int =    file_out.createDataSet("D_stix_int", D_stix_int);
    //h5::DataSet dataset_P_stix_int =    file_out.createDataSet("P_stix_int", P_stix_int);
    //h5::DataSet dataset_R_stix_int =    file_out.createDataSet("R_stix_int", R_stix_int);
    //h5::DataSet dataset_L_stix_int =    file_out.createDataSet("L_stix_int", L_stix_int);
    //h5::DataSet dataset_wmega_e_int =   file_out.createDataSet("wmega_e_int", wmega_e_int);
    //h5::DataSet dataset_mu_ray =        file_out.createDataSet("mu_ray", mu_ray );
    //h5::DataSet dataset_mu_sq =         file_out.createDataSet("mu_sq", mu_sq );
    //h5::DataSet dataset_spsi =          file_out.createDataSet("spsi", spsi );
    //h5::DataSet dataset_cpsi =          file_out.createDataSet("cpsi", cpsi );
    h5::DataSet dataset_kx_ray =        file_out.createDataSet("kx_ray", kx_ray );
    h5::DataSet dataset_kz_ray =        file_out.createDataSet("kz_ray", kz_ray );
    h5::DataSet dataset_kappa_ray =     file_out.createDataSet("kappa_ray", kappa_ray );
    //h5::DataSet dataset_X_stix =        file_out.createDataSet("X_stix", X_stix );
    //h5::DataSet dataset_rho1 =          file_out.createDataSet("rho1", rho1 );
    //h5::DataSet dataset_rho2 =          file_out.createDataSet("rho2", rho2 );
    //h5::DataSet dataset_Byw_sq =        file_out.createDataSet("Byw_sq", Byw_sq );
    //h5::DataSet dataset_fac1 =          file_out.createDataSet("fac1", fac1 );
    //h5::DataSet dataset_Byw =           file_out.createDataSet("Byw", Byw );
    //h5::DataSet dataset_Bxw =           file_out.createDataSet("Bxw", Bxw );
    h5::DataSet dataset_Bzw =           file_out.createDataSet("Bzw", Bzw );
    //h5::DataSet dataset_Exw =           file_out.createDataSet("Exw", Exw );
    //h5::DataSet dataset_Eyw =           file_out.createDataSet("Eyw", Eyw );
    h5::DataSet dataset_Ezw =           file_out.createDataSet("Ezw", Ezw );
    h5::DataSet dataset_Bw_ray =        file_out.createDataSet("Bw_ray", Bw_ray );
    h5::DataSet dataset_w1 =            file_out.createDataSet("w1", w1 );
    h5::DataSet dataset_w2 =            file_out.createDataSet("w2", w2 );
    h5::DataSet dataset_R1 =            file_out.createDataSet("R1", R1 );
    h5::DataSet dataset_R2 =            file_out.createDataSet("R2", R2 );
 

//----------------------------------------------- WRITE HDF5 FILE: DONE -------------------------------------------- //

    auto stop4= std::chrono::high_resolution_clock::now();
    auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
    real time4 = duration4.count()*pow(10,-6);
    std::cout << "\nWriting in h5 file. Time elapsed: " << time4 <<" seconds\n" ;

    return 0;

}

#include "headers/struct_Particles.h"

void Particles::initialize(real eta0, real aeq0, real alpha0, real lamda0, real Ekev0, real Blam0, real zeta0, real deta_dt0, real time0)
{
	//Find momentum from energy.
	real Ejoule0=1.602176487E-16*Ekev0; //Kev to Joule
	real gama0=(Ejoule0/(Constants::m_e*pow(Constants::c,2))) + 1;
	real speed0=sqrt( 1 - (1/pow(gama0,2)) ) * Constants::c;
	//std::cout<<"\nBouncing period estimation: "<< (4*Constants::L_shell*Constants::Re/speed0)*(1.3 - 0.5*sin(aeq0)); //[Orlova1,Shprits2,2011]
	real upar0=speed0*cos(alpha0);
	real uper0=speed0*sin(alpha0);
	real pper0=gama0*Constants::m_e*uper0;  
	real ppar0=gama0*Constants::m_e*upar0;
	real mu_adiabatic_0 = (pper0*pper0) / (2*Constants::m_e*Blam0);
	real Ekin0 = ((gama0-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15; //Joule back to Kev

	//Assign initial state.
	this->lamda.push_back(lamda0);	//PUSH-BACK in RADS
	this->aeq.push_back(aeq0);
	this->eta.push_back(eta0);
	this->alpha.push_back(alpha0);
	this->upar.push_back(upar0);
	this->uper.push_back(uper0);
	this->pper.push_back(pper0);
	this->ppar.push_back(ppar0);
	this->zeta.push_back(zeta0);
	this->time.push_back(time0);
	this->deta_dt.push_back(deta_dt0);
	this->M_adiabatic.push_back(mu_adiabatic_0);	
	this->Ekin.push_back(Ekin0);	

	//Else Invalid
	/*-1<salpha || salpha>1 => domain error
	aeq0=0||180 => salpha0 = 0 => alpha0 = 0 => pper0 = 0 => PROBLEM IN Bell_param, a2(division with 0).
	sin(pi) is actually not zero in C++... */

}

//Member function to push_back new state
void Particles::save_state(int id, real new_lamda, real new_alpha, real new_aeq, real new_ppar, real new_pper, real new_time)
{	//define size
	this->id.push_back(id);
	this->lamda.push_back(new_lamda);      				
	//this->zeta.push_back(new_zeta);
	//this->upar.push_back(new_upar);
	//this->uper.push_back(new_uper);	     
	this->ppar.push_back(new_ppar);		 
	this->pper.push_back(new_pper);
	this->alpha.push_back(new_alpha);	
	this->aeq.push_back(new_aeq);
	//this->deta_dt.push_back(new_deta_dt);
	//this->eta.push_back(new_eta);			
	//this->M_adiabatic.push_back(new_M_adiabatic); 
	//this->Ekin.push_back(new_Ekin); 
	this->time.push_back(new_time);
}	

#include"headers/struct_Telescope.h"

//Constructor. Initialize position.
Telescope::Telescope(real lat, real L_parameter)
{
	L_shell   = L_parameter;
	latitude  = lat;
}
			
bool Telescope::crossing(real p1_lamda, real p2_lamda, real p_L_shell)
{
	bool crossed = false;

	if(p_L_shell == L_shell)
	{
		if( (latitude-p1_lamda)*(latitude-p2_lamda) < 0 ) 
		{ //e.g if particle was below satellite and then above, product would be negative.
			crossed = true;
		}
	}

	return crossed;
}

void Telescope::store(int id, real lamda,real alpha, real time)
{
	this->id.push_back(id);
	this->lamda.push_back(lamda);      				
	//this->upar.push_back(upar);
	//this->uper.push_back(uper);	     
	this->alpha.push_back(alpha);	
	//this->aeq.push_back(aeq);
	//this->eta.push_back(eta);			
	this->time.push_back(time);
}

//Stix tuple, needs wps,wc and returns SDPRL
std::tuple<real, real, real, real, real> stix_parameters(real wce, real wcO, real wcH, real wcHe, real wpse, real wpsO, real wpsH, real wpsHe){ 
//B0mag magnitude of the ambient magnetic field
    real R=1-(wpse/(Constants::w_wave*(Constants::w_wave+wce)))-(wpsH/(Constants::w_wave*(Constants::w_wave+wcH)))-(wpsHe/(Constants::w_wave*(Constants::w_wave+wcHe)))-(wpsO/(Constants::w_wave*(Constants::w_wave+wcO)));
    real L=1-(wpse/(Constants::w_wave*(Constants::w_wave-wce)))-(wpsH/(Constants::w_wave*(Constants::w_wave-wcH)))-(wpsHe/(Constants::w_wave*(Constants::w_wave-wcHe)))-(wpsO/(Constants::w_wave*(Constants::w_wave-wcO)));
    real P=1-(wpse/pow(Constants::w_wave,2))-(wpsH/pow(Constants::w_wave,2))-(wpsHe/pow(Constants::w_wave,2))-(wpsO/pow(Constants::w_wave,2));
    real S = (R+L)/2;
    real D = (R-L)/2 ;
    return std::make_tuple(S, D, P, R, L); //try std::vector for performance
}

//Magnetic dipole field 
real Bmag_dipole(real lamda)                                      //Components of B [Walt,1994]
{   //----calculate the dipole magnetic field strength            //B = sqrt(Br^2 + Bl^2) = B0*(Re/r^3)*slat_term
    //lamda geomagnetic latitude                                  //Where r = Re*L*(clat^2)
    real slat=sin(lamda);                                         //Substituting in Bmag
    real clat=cos(lamda);                   
    real slat_term = sqrt(1 + 3*slat*slat);             
    real Bmag=((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/pow(clat,6);
    return Bmag;  
}

//Solve dispersion relation
std::tuple<real, real, real, real> dispersion(real S, real P, real R, real L, real D){
    //----- solve dispersion relation (find refractive index and k vector)
    //th wave normal angle
    real A=S*sin(Constants::theta0)*sin(Constants::theta0)+P*cos(Constants::theta0)*cos(Constants::theta0);
    real B=R*L*sin(Constants::theta0)*sin(Constants::theta0)+S*P*(1+cos(Constants::theta0)*cos(Constants::theta0));
    real C=P*R*L;                                   //using Lorentz & Maxwell equations with Stix parameters //produced dispersion relation  A*(mu^4) − B*(mu^2) + C = 0
    real mu_sq, mu;
    //std::cout<<"\nB "<< B << " A " << A << " C " << C <<"\nsqrtof "<<(B*B-4*A*C);
    //mu squared can become negative. Probably it has to do with the ion densities. Domain error.
    if(B>=0)  mu_sq=(B-sqrt(B*B-4*A*C))/(2*A);   
    else      mu_sq=(2*C)/(B+sqrt(B*B-4*A*C));      //Conditional forms are only for computational accuracy[Kimura, 1966].
    //mu_sq=(B-sqrt(B*B-4*A*C))/(2*A);
    //std::cout<<"\nmu_sq "<< mu_sq<<" mu "<< mu;
    mu=sqrt(mu_sq);
    real kappa=mu*Constants::w_wave/Constants::c;
    real kx=kappa*sin(Constants::theta0);
    real kz=kappa*cos(Constants::theta0);
    return std::make_tuple(mu, kappa, kx, kz); //try std::vector for performance
}

//Estimate dwh_ds
real dwh_dsf(real w_h, real lamda){         //[Tao et al, 2012]
    real slat = sin(lamda);
    real clat = cos(lamda);
    real slat_term = sqrt(1 + 3*slat*slat);
    real dwh_ds = (3.0*w_h/(Constants::L_shell*Constants::Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat));
    // dwh_ds = (3.0*w_h/(L_shell*Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term))
    return dwh_ds;
}
    
//Estimate Bell parameters //In[10]:                //Bell[1984]
void Bell_params(const real ppar, const real pper, const real Bxw, const real Byw, const real Exw, const real Eyw, const real Ezw, const real kz, const real kx, const real w_h, real gama, real &w1, real &w2, real &wtau_sq, real &R1, real &R2, real &beta) {
    w1=(Constants::q_e/(2*Constants::m_e))*(Bxw+Byw); 
    w2=(Constants::q_e/(2*Constants::m_e))*(Bxw-Byw);
    real wtau0_sq=(w1*kz*pper)/(gama*Constants::m_e);  //Borntik thesis 2.25d
    beta=(kx*pper)/(Constants::m_e*gama*w_h);          //Borntik thesis 2.25a
    real a1=w2/w1;                          //Borntik thesis 2.25f
    real a2=(Constants::q_e*Ezw)/(w1*pper);            //Borntik thesis 2.25g 
    R1=(Exw+Eyw)/(Bxw+Byw);                 //Borntik thesis 2.25h
    R2=(Exw-Eyw)/(Bxw-Byw);                 //Borntik thesis 2.25h
    wtau_sq = (pow((-1),(Constants::m_res-1)) * wtau0_sq * 
            ( jn( ((Constants::m_res)-1), beta ) - 
                a1*jn( ((Constants::m_res)+1) , beta ) +
                gama*a2*jn( Constants::m_res , beta ) )) ;
//Borntik thesis 2.25c

//Bessel function may introduce some inaccuracy(13th decimal)
//All arguments(gama,beta,a1,a2,wtau0_sq) are exactly the same, jn() is not.  
    
//Ji first Kind Bessel function
//integer Order i=(m_res-1)||(m_res+1)||mres
//cylinder function because we have integer order, solving problem in cylindrical coordinate system
//https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_j
//compile with c++17, bessel is defined in <cmath>
//which bessel should i use?
   return; 
}

//Estimate resonant velocity
void vres_f(const real kz, const real w_h, const real alpha, real &v_para_res, real &E_res) {
    real t1 = Constants::w_wave*Constants::w_wave*kz*kz;
    real t2 = pow((Constants::m_res)*w_h, 2)-(Constants::w_wave*Constants::w_wave);
    real t3 = kz*kz + pow((Constants::m_res*w_h),2)/(pow((Constants::c)*cos(alpha),2));
    real direction;                 //positive direction -> counter-streaming particle, negative-> co-streaming
    if((Constants::m_res)==0){
        direction=-1.0*copysign(1,kz);}                             //unsigned m_res                     
    else{                                                           
        direction = copysign(1,kz)*copysign(1,Constants::m_res);}   //signed m_res
    v_para_res = ( direction*sqrt(t1 + (t2*t3)) - Constants::w_wave*kz) / t3;
    real v_tot_res = v_para_res / cos(alpha);
    E_res = Constants::e_el*(1.0/sqrt( 1-(v_tot_res*v_tot_res/(Constants::c*Constants::c)) ) - 1 );  //not sure about that... taken from Bortnik code
    return;
}

//Compute field components.([Bell 1984] 𝑧̂ ||𝐵0→ and 𝑥̂  pointing towards higher L-shells)
void whistlers(int64_t p, int64_t i, real mu, real P, real D, real S, real kz, real &Bxw,real &Byw,real &Bzw, real &Exw,real &Eyw,real &Ezw)
{

    real mu_sq = pow(mu,2); //mu*mu inaccurate
    real theta = Constants::theta0;
    real fac1 = P-(mu_sq*pow(sin(theta),2)) ;
    Byw = Constants::By_wave;

    Bxw = std::abs((-(D*fac1)/(P*(S-mu_sq)))*Byw); 
    Bzw = std::abs(((D*sin(theta)*fac1)/(P*cos(theta)*(S-(mu_sq)))*Byw));
    Exw = std::abs((((Constants::c)*fac1)/(mu*P*cos(theta))*Byw));
    Eyw = std::abs(((D*(Constants::c)*fac1)/(mu*P*cos(theta)*(mu_sq-S)))*Byw);
    Ezw = std::abs((-((Constants::c)*mu*sin(theta))/P)*Byw);

    //real Phi  = Constants::w_wave*time-kz*zeta; 
    //Exw = -Exw.at(i)*sin(Phi);  
    //Eyw =  Eyw.at(i)*cos(Phi);
    //Ezw = -Ezw.at(i)*sin(Phi);
    //Bxw =  Bxw.at(i)*cos(Phi);  
    //Byw =  Byw.at(i)*sin(Phi);
    //Bzw = -Bzw.at(i)*cos(Phi);

    return;
}






//For Ray Tracing


//Return Factors, aeq of runge kutta and kz vector, only when particle in packet.
void f_packet (real &Fpar, real &Fper, real &Ftheta, real &aeq_rk, real &kz, const real pper_tmp, const real ppar_tmp, const real eta_tmp, const real aeqsu_tmp, const real alpha_tmp, const real gama, const real w_h, const real p_mag, const real kx_ray, const real kz_ray, const real kappa_ray, const real Bw_ray, const real Bzw, const real Ezw, const real w1, const real w2, const real R1, const real R2 )
{
    //Takes and returns same value. Change that
    kz       =    kz_ray; 
    
    //Using wave's values where the particle is confined...
    

    //std::cout<<"\npper_tmp "<<pper_tmp<<"\nppar_tmp "<<ppar_tmp<<"\neta_tmp "<<eta_tmp<<"\naeqsu_tmp "<<aeqsu_tmp<<"\nalpha_tmp "<<alpha_tmp<<"\ngama "<<gama<<"\nw_h "<<w_h<<"\np_mag "<<p_mag<<"\nkx_ray "<<kx_ray<<"\nkz_ray "<<kz_ray<<"\nkappa_ray"<<kappa_ray<<"\nBw_ray "<<Bw_ray<<"\nEzw "<<Ezw<<"\nw1 "<<w1<<"\nw2 "<<w2<<"\nR1 "<<R1<<"\nR2 "<<R2<<"\n";


    //Calculate beta 
    real beta = (kx_ray*pper_tmp)/(Constants::m_e*gama*w_h); //[Bortnic thesis]. No need to return.
    //std::cout<<"\nbeta "<<beta;     //It's negative-> domain error for std::cyl_bessel_j().

    
    //Factor calculations for ppar ,pper and theta
    real EL = R1*((2*Constants::m_e)/Constants::q_e)*w1;
    real ER = R2*((2*Constants::m_e)/Constants::q_e)*w2;
    real BR = ((2*Constants::m_e)/Constants::q_e)*w1;
    real BL = ((2*Constants::m_e)/Constants::q_e)*w2;
    
    Fpar = -(pow((-1),(Constants::m_res-1))*(-Constants::q_e*(Ezw*jn( (Constants::m_res), beta)+(pper_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res-1), beta)
                   -(pper_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res+1), beta))));
                   
    Fper = (-Constants::q_e*(ER*jn( (Constants::m_res-1), beta)+EL*jn( (Constants::m_res+1), beta)
                  -(ppar_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res-1), beta)
                  +(ppar_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res+1), beta)));

    Ftheta = (-Constants::q_e*(ER*jn( (Constants::m_res-1), beta) - EL*jn( (Constants::m_res+1), beta) 
                  -(ppar_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res-1), beta) 
                  -(ppar_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res+1), beta)
                  +(pper_tmp/(gama*Constants::m_e))*Bzw*jn( (Constants::m_res), beta)));
    
    //Calculate Equatorial P.A change.
    aeq_rk =((Constants::q_e*Bw_ray)/(pow(p_mag,2)))*(tan(aeqsu_tmp)/tan(alpha_tmp))*(((Constants::w_wave/kappa_ray)-(ppar_tmp/(gama*Constants::m_e)))*ppar_tmp-(pow(pper_tmp,2)/(gama*Constants::m_e)))*sin(eta_tmp);
    return;
}

//Returns quantities needed, regardless WPI.
void f_always(real &p_mag, real &gama, real &w_h, real &dwh_ds, const real lamda_tmp, const real ppar_tmp, const real pper_tmp)
{   
    
    p_mag  =  sqrt(ppar_tmp*ppar_tmp + pper_tmp*pper_tmp);
    gama   =  sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
    //same with gama = sqrt( 1+(u/c)^2 )
    real slat = sin(lamda_tmp);   //Bmag and dwh_ds needs                                
    real clat = cos(lamda_tmp);                   
    real slat_term = sqrt(1 + 3*slat*slat); 
    
    real Bmag = ((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/(pow(clat,6));
    w_h = ((Constants::q_e*Bmag)/Constants::m_e);
    dwh_ds = (3.0*w_h/(Constants::L_shell*Constants::Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat)); // [Tao et al, 2012]
    return;
}

//Estimate resonant velocity
real vres_f(real w_wave, real kz, real w_h, real alpha) 
{
    real t1 = w_wave*w_wave*kz*kz;
    real t2 = pow((Constants::m_res)*w_h, 2)-(w_wave*w_wave);
    real t3 = kz*kz + pow((Constants::m_res*w_h),2)/(pow((Constants::c)*cos(alpha),2));
    real direction;                                                 //Positive direction when counter-streaming particle
    if((Constants::m_res)==0)                                       //Negative is co-streaming
    { direction=-1.0*copysign(1,kz); }                              //Unsigned m_res                     
    else                                                            
    { direction = copysign(1,kz)*copysign(1,Constants::m_res); }    //Signed m_res
    real v_para_res = ( direction*sqrt(t1 + (t2*t3)) - w_wave*kz) / t3;
    //real v_tot_res = v_para_res / cos(alpha);
    //real E_res = Constants::e_el*(1.0/sqrt( 1-(v_tot_res*v_tot_res/(Constants::c*Constants::c)) ) - 1 ); //not used
    return v_para_res; //Return by reference, global variables, temp.
}

#include "headers/li_wpi.h"

#include "headers/is_in_packet.h"

//Check if the particle is within the wave's latitude.
//Return wave's index where there's matching latitude(WPI). Will be used for wave's data.
int is_in_packet(real min_lat, real max_lat, real lamda_tmp, int i, std::vector<real> &wave_lat)
{
	std::vector <real>::iterator ptr, requested_index;    //To find requested index
	int index;
	real diff, min_diff;                                  //To find minimum difference between latitudes

          
    if(min_lat < (lamda_tmp*Constants::R2D) && (lamda_tmp*Constants::R2D) < max_lat)
    {       
		min_diff = 100 ; 
		for(ptr = wave_lat.begin() + i; ptr < wave_lat.begin() + (Constants::puls_dur + i); ptr++)
        {
       	    diff = std::abs(*ptr - lamda_tmp*Constants::R2D) ;        //Return abs of difference
       	    if(diff < min_diff) 
       	    {          
			    requested_index = ptr; //Current requested index
			   	min_diff = diff;   	   //Current minimum                                                             
       	    }
       	}
	index = distance(wave_lat.begin(), requested_index); //Final requested index
	//std::cout<<"\nindex "<<int(index-Constants::puls_dur)<<"\n";
	}
	else
	{
		index = -1;	//Particle not in same latitude. No WPI.
	}

   
	return index ;
}


void wpi_ray(real p, Particles &single, Telescope &ODPT)
{
//---------------------------------------------------- READ RAY HDF5 ----------------------------------------------------//
    //read_vector() is a function to read HDF5 dataset as vectors. 
    static std::vector <real> lat_int       =   read_vector("lat_int",       "h5files/interpolated_ray.h5");
    static std::vector <real> kx_ray        =   read_vector("kx_ray",        "h5files/interpolated_ray.h5");    
    static std::vector <real> kz_ray        =   read_vector("kz_ray",        "h5files/interpolated_ray.h5");   
    static std::vector <real> kappa_ray     =   read_vector("kappa_ray",     "h5files/interpolated_ray.h5");       
    static std::vector <real> Bzw           =   read_vector("Bzw",           "h5files/interpolated_ray.h5");
    static std::vector <real> Ezw           =   read_vector("Ezw",           "h5files/interpolated_ray.h5");
    static std::vector <real> Bw_ray        =   read_vector("Bw_ray",        "h5files/interpolated_ray.h5");    
    static std::vector <real> w1            =   read_vector("w1",            "h5files/interpolated_ray.h5");
    static std::vector <real> w2            =   read_vector("w2",            "h5files/interpolated_ray.h5");
    static std::vector <real> R1            =   read_vector("R1",            "h5files/interpolated_ray.h5");
    static std::vector <real> R2            =   read_vector("R2",            "h5files/interpolated_ray.h5");
//---------------------------------------------------- ASSIGN OBJECT VALUES ----------------------------------------------------//
                                           //Erase first element which was the last state of the previous noWPI simulation.
    real lamda    =  single.lamda.at(0);    single.lamda.erase(single.lamda.begin());
    real ppar     =  single.ppar.at(0);     single.ppar.erase(single.ppar.begin()); 
    real pper     =  single.pper.at(0);     single.pper.erase(single.pper.begin()); 
    real eta      =  single.eta.at(0);      single.eta.erase(single.eta.begin()); 
    real alpha    =  single.alpha.at(0);    single.alpha.erase(single.alpha.begin()); 
    real aeq      =  single.aeq.at(0);      single.aeq.erase(single.aeq.begin());
    real time     =  single.time.at(0);     single.time.erase(single.time.begin());
    //real zeta     =  single.zeta.at(0);           single.zeta.erase(single.zeta.begin()); 
    //real upar     =  single.upar.at(0);           single.upar.erase(single.upar.begin()); 
    //real uper     =  single.uper.at(0);           single.uper.erase(single.uper.begin());
    //real deta_dt  =  single.deta_dt.at(0);        single.deta_dt.erase(single.deta_dt.begin());
    //real Ekin     =  single.Ekin.at(0);           single.Ekin.erase(single.Ekin.begin());
    //real M_adiabatic = single.M_adiabatic.at(0);  single.M_adiabatic.erase(single.M_adiabatic.begin());
//------------------------------------------------- LOOP DECLARATIONS -------------------------------------------------//
    int index;                              //To find minimum difference between latitudes
    int i=0;

    real min_lat,max_lat;
    real p_mag,gama,w_h,dwh_ds,kz,Fpar,Fper,Ftheta;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
    real new_lamda, new_aeq, new_ppar;
    bool trapped = 1;                       //Particles trapped in Earth's magnetic field.

    //std::cout.precision(64);                //Output 16 decimal precise
	//std::cout<<std::scientific;		        //For e notation representation
//----------------------------------------------------- WPI -----------------------------------------------------------//

    while(i<Constants::Nsteps_wpi)          
    {   
        //Take the data from interpolation and return (pulse_dur) of them, moved by "i" each iteration.
        min_lat=*min_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Minimum lat of wave(in pulse duration).
        max_lat=*max_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Max lat of wave.
        //std::cout<<"\n" << lamda*Constants::R2D << " " << min_lat << "" << max_lat << " "<< lat_int.at(i) << " " << lat_int.back(); 
        
        f_always(p_mag, gama, w_h, dwh_ds, lamda, ppar, pper); 
        kz = Ftheta = Fpar = Fper = q1 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda, i, lat_int);  //is_in_packet returns "if and where" WPI happens. 
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q1, kz, pper, ppar, eta, aeq, alpha, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes(k1,  l1,  m1,  n1,  o1,  p1, ppar, pper, alpha, lamda, eta, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz);
        //std::cout <<"\np_mag "<<p_mag<<"\ndwh_ds "<< dwh_ds<<"\nkz "<< kz  <<"\nFpar "<< Fpar <<"\nFper "<< Fper <<"\nFtheta "<< Ftheta <<"\n";  
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn1 " << n1<< "\no1 " << o1 << "\np1 " << p1 <<"\nq1 " << q1 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, lamda+0.5*Constants::h*o1, ppar+0.5*Constants::h*l1, pper+0.5*Constants::h*m1);
        kz = Ftheta = Fpar = Fper = q2 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda+0.5*Constants::h*o1, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q2, kz, pper+0.5*Constants::h*m1, ppar+0.5*Constants::h*l1, eta+0.5*Constants::h*n1, aeq+0.5*Constants::h*q1, alpha+0.5*Constants::h*p1, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k2, l2, m2, n2, o2, p2, ppar + 0.5*Constants::h*l1, pper + 0.5*Constants::h*m1, alpha + 0.5*Constants::h*p1, lamda + 0.5*Constants::h*o1, eta + 0.5*Constants::h*n1, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz ); 
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\nq2 "<< q2 <<"\n";;
        
        f_always(p_mag, gama, w_h, dwh_ds, lamda+0.5*Constants::h*o2, ppar+0.5*Constants::h*l2, pper+0.5*Constants::h*m2);
        kz = Ftheta = Fpar = Fper = q3 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda+0.5*Constants::h*o2, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q3, kz, pper+0.5*Constants::h*m2, ppar+0.5*Constants::h*l2, eta+0.5*Constants::h*n2, aeq+0.5*Constants::h*q2, alpha+0.5*Constants::h*p2, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k3, l3, m3, n3, o3, p3, ppar + 0.5*Constants::h*l2, pper + 0.5*Constants::h*m2, alpha + 0.5*Constants::h*p2, lamda + 0.5*Constants::h*o2, eta + 0.5*Constants::h*n2, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz );   
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\nq3 "<< q3 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, lamda+Constants::h*o3, ppar+Constants::h*l3, pper+Constants::h*m3);
        kz = Ftheta = Fpar = Fper = q4 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda+Constants::h*o3, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q4, kz, pper+Constants::h*m3, ppar+Constants::h*l3, eta+Constants::h*n3, aeq+Constants::h*q3, alpha+Constants::h*p3, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k4, l4, m4, n4, o4, p4, ppar + Constants::h*l3, pper + Constants::h*m3, alpha + Constants::h*p3, lamda + Constants::h*o3, eta + Constants::h*n3, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz ); 
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn4 " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\nq4 " << q4 <<"\n";           


        //Check validity:
        new_lamda = lamda + (Constants::h/6)*(o1+2*o2+2*o3+o4); //Approximate new lamda first
        if(std::isnan(new_lamda)) { std::cout<<"\nParticle "<<p<<" breaks"; break; }
        //Check crossing:
        #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
        {                    //This slows down the parallel process, introduces bad scalling 8+ cores. Detecting first and storing in the end demands more memory per process.
            if( ODPT.crossing(new_lamda*Constants::R2D, lamda*Constants::R2D, Constants::L_shell) )	 
            {	//Check crossing.								
                //std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
                ODPT.store( p, lamda, alpha, time); //Store its state(it's before crossing the satellite!).		        	
            }
        }
        //Check precipitation:
        new_aeq = aeq + (Constants::h/6)*(q1+2*q2+2*q3+q4);
        if(new_aeq<Constants::alpha_lc) //If particle's equator P.A is less than the loss cone angle for this L_shell, then particle is not trapped. Minimum allowable trapping altitude 100km.
        {
            trapped = 0;
        }
        //If it's not trapped and it's about to bounce --> Precipitation
        new_ppar = ppar + (Constants::h/6)*(l1+2*l2+2*l3+l4);
        if(!trapped && (ppar*new_ppar<0) ) //Would bounce if ppar is about to change sign
        {   
            //To save states of precipitating particles:
            #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
            {
                single.save_state(p, lamda, alpha, aeq, ppar, pper, time);
                std::cout<<"\n\nParticle "<<p<<" escaped with ppar "<<ppar<< " new_ppar would be "<<new_ppar<<" pper " << pper<< " eta " << eta << " lamda " <<lamda<< " alpha "<< alpha << " aeq " <<aeq<< " at time " << time ;
            }
            break;
        }
        //Next step:
        new_values_RK4(lamda, ppar, pper, eta, alpha, aeq, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, o1, o2, o3, o4, p1, p2, p3, p4, q1, q2, q3, q4);
        time  = time + Constants::h; 
        i++;  
        //std::cout<<"\n\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda<< "\nalpha "<< alpha << "\naeq " <<aeq ;
    }


}

#include "headers/bell_nowpi.h"

//Adiabatic motion.
void nowpi(int p, Particles &single, Telescope &ODPT)
{
    real lamda   =  single.lamda.at(0);
    //real zeta    =  single.zeta.at(0); 
    real ppar    =  single.ppar.at(0); 
    real pper    =  single.pper.at(0); 
    real alpha   =  single.alpha.at(0); 
    real aeq     =  single.aeq.at(0); 
    //real upar    =  single.upar.at(0); 
    //real uper    =  single.uper.at(0);
    //real Ekin    =  single.Ekin.at(0);
    real time    =  single.time.at(0);

    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamda;
    real w_h, dwh_ds, Bmag, p_mag, gama;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,o1,o2,o3,o4,p1,p2,p3,p4;

    //Objects for each specie.
    Species electron(Constants::m_e,  Constants::q_e, 1); 
    Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
    Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
    Species helium  (Constants::m_He, Constants::q_i, 0.054);
    
	//std::cout.precision(64);			//Output 16 decimal precise
	//std::cout<<std::scientific;		//For e notation representation

    int i=0;
    
    while(i<Constants::Nsteps_nowpi) 
    {

        Bmag=Bmag_dipole(lamda);   
        w_h = electron.wc(Bmag); //Cyclotron frequency.
        dwh_ds=dwh_dsf(w_h,lamda);
        p_mag = sqrt(ppar*ppar+pper*pper);
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-1//#################################################################################################################################################################################################################################################################################################
        slopes(k1, l1, m1, o1, p1, ppar, pper, lamda, w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 <<"\no1 " << o1 << "\np1 " << p1 <<"\n";	
        

        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o1);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o1);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l1)*(ppar+0.5*(Constants::h)*l1)+(pper+0.5*(Constants::h)*m1)*(pper+0.5*(Constants::h)*m1));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-2//#################################################################################################################################################################################################################################################################################################
        slopes(k2, l2, m2, o2, p2, ppar+(0.5*l1*Constants::h), pper+(0.5*m1*Constants::h), lamda+(0.5*o1*Constants::h),w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 <<"\no2 " << o2 << "\np2 " << p2 << "\n";
        

        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o2);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o2);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l2)*(ppar+0.5*(Constants::h)*l2)+(pper+0.5*(Constants::h)*m2)*(pper+0.5*(Constants::h)*m2));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-3//#################################################################################################################################################################################################################################################################################################
        slopes(k3, l3, m3, o3, p3, ppar+(0.5*l2*Constants::h), pper+(0.5*m2*Constants::h), lamda+(0.5*o2*Constants::h),w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3  << "\no3 " << o3 << "\np3 " << p3 << "\n";


        Bmag=Bmag_dipole(lamda+(Constants::h)*o3);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+(Constants::h)*o3);
        p_mag = sqrt((ppar+(Constants::h)*l3)*(ppar+(Constants::h)*l3)+(pper+(Constants::h)*m3)*(pper+(Constants::h)*m3));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-4//#################################################################################################################################################################################################################################################################################################																								
        slopes(k4, l4, m4, o4, p4, ppar+(l3*Constants::h), pper+(m3*Constants::h), lamda+(o3*Constants::h), w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 <<"\no4 " << o4 << "\np4 " << p4 << "\n";


        //Approximate new lamda first, to check if particle crosses satellite.
        new_lamda = lamda + ((Constants::h)/6)*(o1+2*o2+2*o3+o4);
        #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
        {                    //This slows down the parallel process, introduces bad scalling 8+ cores. Detecting first and storing in the end demands more memory per process.
            if( ODPT.crossing(new_lamda*Constants::R2D, lamda*Constants::R2D, Constants::L_shell) )	 
            {	//Check crossing.								
                //std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
                ODPT.store( p, lamda, alpha, time); //Store its state(it's before crossing the satellite!).		        	
            }
        }

        //Next step:
        new_values_RK4(lamda, ppar, pper, alpha, l1, l2, l3, l4, m1, m2, m3, m4, o1, o2, o3, o4, p1, p2, p3, p4);
        time  = time + Constants::h; 
        i++;  

		//To save states:
		//single.save_state(aeq,alpha,lamda,deta_dt,time);
        //std::cout<<"\n\nalpha "<<alpha << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda<< "\naeq " <<aeq ;

        //Stop at equator:
        //if(eql_dstr[p].lamda.at(i)>0) {	
        //	break;}	
    }
    

    //Save last state to continue the simulation with wave if needed. 
    //single.save_state(lamda,alpha,aeq,ppar,pper,time);
    single.lamda.front() = lamda;
    single.alpha.front() = alpha;
    single.aeq.front()   = aeq;
    single.ppar.front()  = ppar;
    single.pper.front()  = pper;
    single.eta.front()   = Constants::eta0;
    single.time.front()  = time;
    
}


#include "headers/bell_wpi.h"

//For wave-particle interaction
void wpi(int p, Particles &single, Telescope &ODPT)
{
    
	std::cout.precision(64);			//Output 16 decimal precise
	std::cout<<std::scientific;		    //For e notation representation
    
    real lamda    =  single.lamda.at(0);
    //real zeta     =  single.zeta.at(0); 
    real ppar     =  single.ppar.at(0); 
    real pper     =  single.pper.at(0); 
    real alpha    =  single.alpha.at(0); 
    real aeq      =  single.aeq.at(0); 
    //real upar     =  single.upar.at(0); 
    //real uper     =  single.uper.at(0);
    //real deta_dt  =  single.deta_dt.at(0);
    //real M_adiabatic = single.M_adiabatic.at(0);
    //real Ekin     =  single.Ekin.at(0);
    real eta      =  single.eta.at(0); 
    real time     =  single.time.at(0);
    std::cout<<"\n\nParticle "<<p<<" at alpha "<<alpha << "\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda<< "\naeq " <<aeq ;

    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamda = lamda;
    real ns_e,w_h, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He;
    real Bmag;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
    real l1_old,l2_old,l3_old,l4_old,m1_old,m2_old,m3_old,m4_old,n1_old,n2_old,n3_old,n4_old,o1_old,o2_old,o3_old,o4_old,p1_old,p2_old,p3_old,p4_old,q1_old,q2_old,q3_old,q4_old;
    real gama,w1,w2,R1,R2,beta,wtau_sq;
    real S,D,P,R,L,mu,dwh_ds,kappa,kx,kz;
    //real vres, Eres;
    real Bxwc,Bzwc,Bywc,Exwc,Eywc,Ezwc,Bwc;
    real p_mag;
    //Tuples for WPI
    std::tuple<real, real, real, real, real> stix;
    std::tuple<real, real, real, real> disp;
    //Objects for each specie.
    Species electron(Constants::m_e,  Constants::q_e, 1); 
    Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
    Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
    Species helium  (Constants::m_He, Constants::q_i, 0.054);
    

    int i=0;

    while(i<Constants::Nsteps_wpi) 
    {
        Bmag=Bmag_dipole(lamda);
        ns_e = electron.density(lamda); ns_O = oxygen.density(lamda);  ns_H = hydrogen.density(lamda); ns_He = helium.density(lamda);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);  wc_He = helium.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda);
        p_mag = sqrt(ppar*ppar+pper*pper);
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix); //<get> efficiency comperable to accessing a member of a struct
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar,pper,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        //vres_f(kz,w_h,alpha,vresz,Eres); //Called only once in first step...
        if(std::isnan(mu)) //mu becomes nan first
        {
            if(i==0) 
            { 
                std::cout<<"\nParticle "<<p<<" broke at first step."<<std::endl; break;
            }
            else
            {
                //Move to next step to avoid nan. Use old slopes(previous step).
                std::cout<<"\nParticle "<<p<<" at step: "<< i <<" continues."<<std::endl;
                new_values_RK4(lamda, ppar, pper, eta, alpha, aeq, l1_old, l2_old, l3_old, l4_old, m1_old, m2_old, m3_old, m4_old, n1_old, n2_old, n3_old, n4_old, o1_old, o2_old, o3_old, o4_old, p1_old, p2_old, p3_old, p4_old, q1_old, q2_old, q3_old, q4_old);
                time  = time + Constants::h; 
                i++;
                continue;  
            }
        }
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k1, l1, m1, n1, o1, p1, q1, ppar, pper, lamda, eta, alpha, aeq, p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\nq1 " << q1 <<"\n";	


        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o1);
        ns_e = electron.density(lamda+0.5*(Constants::h)*o1); ns_O = oxygen.density(lamda+0.5*(Constants::h)*o1);  ns_H = hydrogen.density(lamda+0.5*(Constants::h)*o1); ns_He = helium.density(lamda+0.5*(Constants::h)*o1);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o1);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l1)*(ppar+0.5*(Constants::h)*l1)+(pper+0.5*(Constants::h)*m1)*(pper+0.5*(Constants::h)*m1));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+0.5*(Constants::h)*l1,pper+0.5*(Constants::h)*m1,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        if(std::isnan(mu)) //mu becomes nan first
        {
            if(i==0) 
            { 
                std::cout<<"\nParticle "<<p<<" broke at first step."<<std::endl; break;
            }
            else
            {
                //Move to next step to avoid nan. Use old slopes(previous step).
                std::cout<<"\nParticle "<<p<<" at step: "<< i <<" continues."<<std::endl;
                new_values_RK4(lamda, ppar, pper, eta, alpha, aeq, l1_old, l2_old, l3_old, l4_old, m1_old, m2_old, m3_old, m4_old, n1_old, n2_old, n3_old, n4_old, o1_old, o2_old, o3_old, o4_old, p1_old, p2_old, p3_old, p4_old, q1_old, q2_old, q3_old, q4_old);
                time  = time + Constants::h; 
                i++;
                continue;  
            }
        }
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k2, l2, m2, n2, o2, p2, q2, ppar+(0.5*l1*Constants::h), pper+(0.5*m1*Constants::h), lamda+(0.5*o1*Constants::h), eta+(0.5*n1*Constants::h), alpha+(0.5*p1*Constants::h), aeq+(0.5*q1*Constants::h), p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 <<"\nq2 "<< q2 <<"\n";
        

        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o2);
        ns_e = electron.density(lamda+0.5*(Constants::h)*o2); ns_O = oxygen.density(lamda+0.5*(Constants::h)*o2);  ns_H = hydrogen.density(lamda+0.5*(Constants::h)*o2); ns_He = helium.density(lamda+0.5*(Constants::h)*o2);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o2);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l2)*(ppar+0.5*(Constants::h)*l2)+(pper+0.5*(Constants::h)*m2)*(pper+0.5*(Constants::h)*m2));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+0.5*(Constants::h)*l2,pper+0.5*(Constants::h)*m2,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);  
        if(std::isnan(mu)) //mu becomes nan first
        {
            if(i==0) 
            { 
                std::cout<<"\nParticle "<<p<<" broke at first step."<<std::endl; break;
            }
            else
            {
                //Move to next step to avoid nan. Use old slopes(previous step).
                std::cout<<"\nParticle "<<p<<" at step: "<< i <<" continues."<<std::endl;
                new_values_RK4(lamda, ppar, pper, eta, alpha, aeq, l1_old, l2_old, l3_old, l4_old, m1_old, m2_old, m3_old, m4_old, n1_old, n2_old, n3_old, n4_old, o1_old, o2_old, o3_old, o4_old, p1_old, p2_old, p3_old, p4_old, q1_old, q2_old, q3_old, q4_old);
                time  = time + Constants::h; 
                i++;
                continue;  
            }
        }
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k3, l3, m3, n3, o3, p3, q3, ppar+(0.5*l2*Constants::h), pper+(0.5*m2*Constants::h), lamda+(0.5*o2*Constants::h), eta+(0.5*n2*Constants::h), alpha+(0.5*p2*Constants::h), aeq+(0.5*q2*Constants::h), p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 <<"\nq3 "<< q3 <<"\n";
        

        Bmag=Bmag_dipole(lamda+Constants::h*o3);
        ns_e = electron.density(lamda+Constants::h*o3); ns_O = oxygen.density(lamda+Constants::h*o3);  ns_H = hydrogen.density(lamda+Constants::h*o3); ns_He = helium.density(lamda+Constants::h*o3);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,lamda+Constants::h*o3);
        p_mag = sqrt((ppar+Constants::h*l3)*(ppar+Constants::h*l3)+(pper+Constants::h*m3)*(pper+Constants::h*m3));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+Constants::h*l3,pper+Constants::h*m3,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);  
        if(std::isnan(mu)) //mu becomes nan first
        {
            if(i==0) 
            { 
                std::cout<<"\nParticle "<<p<<" broke at first step."<<std::endl; break;
            }
            else
            {
                //Move to next step to avoid nan. Use old slopes(previous step).
                std::cout<<"\nParticle "<<p<<" at step: "<< i <<" continues."<<std::endl;
                new_values_RK4(lamda, ppar, pper, eta, alpha, aeq, l1_old, l2_old, l3_old, l4_old, m1_old, m2_old, m3_old, m4_old, n1_old, n2_old, n3_old, n4_old, o1_old, o2_old, o3_old, o4_old, p1_old, p2_old, p3_old, p4_old, q1_old, q2_old, q3_old, q4_old);
                time  = time + Constants::h; 
                i++;
                continue;  
            }
        }
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k4, l4, m4, n4, o4, p4, q4, ppar+l3*Constants::h, pper+m3*Constants::h, lamda+o3*Constants::h, eta+n3*Constants::h, alpha+p3*Constants::h, aeq+q3*Constants::h, p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\nq4 "<< q4 <<"\n";
       
       
        //Approximate new lamda first, to check if particle crosses satellite.
        new_lamda = lamda + (Constants::h/6)*(o1+2*o2+2*o3+o4);
        #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
        {                    //This slows down the parallel process, introduces bad scalling 8+ cores. Detecting first and storing in the end demands more memory per process.
            if( ODPT.crossing(new_lamda*Constants::R2D, lamda*Constants::R2D, Constants::L_shell) )	 
            {	//Check crossing.								
                //std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
                ODPT.store( p, lamda, alpha, time); //Store its state(it's before crossing the satellite!).		        	
            }
        }
        
        // Old slope values kept in memory to encounter the NAN case
        // isnan(mu) --> step once more using old slopes to reach to a valid state
        // i.e when NAN, step of h becomes 2*h, 3*h ... until valid.
        l1_old = l1 ; l2_old = l2 ; l3_old = l3 ; l4_old = l4 ;  
        m1_old = m1 ; m2_old = m2 ; m3_old = m3 ; m4_old = m4 ;  
        n1_old = n1 ; n2_old = n2 ; n3_old = n3 ; n4_old = n4 ;  
        o1_old = o1 ; o2_old = o2 ; o3_old = o3 ; o4_old = o4 ;  
        p1_old = p1 ; p2_old = p2 ; p3_old = p3 ; p4_old = p4 ;  
        q1_old = q1 ; q2_old = q2 ; q3_old = q3 ; q4_old = q4 ; 

        //Next step:
        new_values_RK4(lamda, ppar, pper, eta, alpha, aeq, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, o1, o2, o3, o4, p1, p2, p3, p4, q1, q2, q3, q4);
        time  = time + Constants::h; 
        i++;  

		//To save states:
		//single.save_state(aeq,alpha,lamda,deta_dt,time);
        //std::cout<<"\n\nalpha "<<alpha << "\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda<< "\naeq " <<aeq ;

        //Stop at equator:
        //if(eql_dstr[p].lamda.at(i)>0) {	
        //	break;}	
    }

}


int main(int argc, char **argv)
{
	if(argc!=2)
	{
		std::cout<<"Error. Argc should be 2. Set second argv from the list:(bell, li_ray)."<<std::endl;
		return EXIT_FAILURE;
	}
	std::cout<<"\n"<<Constants::alpha_lc*Constants::R2D;
	
	//Position of the Particle Telescope.		
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell);		
	//Single particle struct.
	Particles single; 
	//Vector of structs for particle distribution.
	std::vector<Particles> dstr(Constants::population, single);	

//------------------------------------------------------------READ DISTRIBUTION FROM H5 FILE --------------------------------------------------------------//
	h5::File dstr_file("h5files/distribution_5000.h5", h5::File::ReadOnly);
	//Vectors to save temporarily
	std::vector<real> lamda_0, alpha_0, aeq_0, ppar_0, pper_0, upar_0, uper_0, Ekin_0, time_0, zeta_0, eta_0, deta_dt_0, M_adiabatic_0;
	//Read dataset from h5file.
	h5::DataSet data_lat 	     = dstr_file.getDataSet("lat");
	h5::DataSet data_aeq	     = dstr_file.getDataSet("aeq");
	h5::DataSet data_alpha 		 = dstr_file.getDataSet("alpha");
	h5::DataSet data_upar 		 = dstr_file.getDataSet("upar");
	h5::DataSet data_uper 		 = dstr_file.getDataSet("uper");
	h5::DataSet data_ppar 		 = dstr_file.getDataSet("ppar");
	h5::DataSet data_pper 		 = dstr_file.getDataSet("pper");
	h5::DataSet data_eta 		 = dstr_file.getDataSet("eta");
	h5::DataSet data_zeta 		 = dstr_file.getDataSet("zeta");
	h5::DataSet data_time 		 = dstr_file.getDataSet("time");
	h5::DataSet data_deta_dt     = dstr_file.getDataSet("deta_dt");
	h5::DataSet data_M_adiabatic = dstr_file.getDataSet("M_adiabatic");
	h5::DataSet data_Ekin 		 = dstr_file.getDataSet("Ekin");
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
		dstr[p].lamda.push_back(lamda_0.at(p));
		dstr[p].alpha.push_back(alpha_0.at(p));  
		dstr[p].aeq.push_back(aeq_0.at(p));
		dstr[p].ppar.push_back(ppar_0.at(p));
		dstr[p].pper.push_back(pper_0.at(p));
		dstr[p].upar.push_back(upar_0.at(p));
		dstr[p].uper.push_back(uper_0.at(p));
		dstr[p].Ekin.push_back(Ekin_0.at(p));
		dstr[p].time.push_back(time_0.at(p));
		dstr[p].zeta.push_back(zeta_0.at(p));
		dstr[p].eta.push_back(eta_0.at(p));
		dstr[p].deta_dt.push_back(deta_dt_0.at(p));
		dstr[p].M_adiabatic.push_back(M_adiabatic_0.at(p));
	}
	std::cout<<"\nParticle population: "<< dstr.size()<<std::endl;


//--------------------------------------------------------------------SIMULATION------------------------------------------------------------------------//
	int realthreads;   
	real wtime = omp_get_wtime();
	std::string s1("bell");
	std::string s2("li_ray");

	//---NOWPI---//
	if(Constants::t_nowpi!=0)
	{
		std::cout<<"\n\n"<<Constants::t_nowpi<<" sec NoWPI Simulation using Bell formulas"<<std::endl;
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
					//std::cout<<"\nBouncing particle "<<p<<" "<<id<<std::flush;
					//Void Function for particle's motion. Involves RK4 for Nsteps. 
					//Detected particles are saved in ODPT object, which is passed here by reference.
					nowpi(p, dstr[p], ODPT);
					//Inside nowpi -> Last states become firsts to continue the simulation 
					//Then wpi
				}
			
		}	
		std::cout<<"\n"<<"Joined"<<std::endl;
		real time1 = omp_get_wtime()-wtime;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time1<<std::endl;
	}


	//---WPI---//
	if(Constants::t_wpi!=0)
	{
		if(!(s1.compare(argv[1])))
		{
			std::cout<<"\n\n"<<Constants::t_wpi<<" sec WPI Simulation using Bell formulas. Wave magnitude(T): "<<Constants::By_wave<<std::endl;
			std::cout<<"Execution time estimation for 8 THREAD run: "<<(Constants::population*0.036/60) * Constants::t_wpi <<" minutes."<<std::endl;
		}
		else if(!(s2.compare(argv[1])))
		{
			std::cout<<"\n\n"<<Constants::t_wpi<<" sec ray tracing WPI Simulation using Li formulas."<<std::endl;
			std::cout<<"Execution time estimation for 20 THREAD run: "<<(Constants::population*0.1/60) * Constants::t_wpi <<" minutes."<<std::endl;
		}
		else
		{ 
			std::cout<<"\n\nArgument variable doesn't match any of the program's possible implementations.\nTry nowpi, wpi, or li_ray as the second argument variable.\n"<<std::endl;
			return EXIT_FAILURE;	
		}

		std::cout<<"\nForked..."<<std::endl;
		//---PARALLELISM Work sharing---//
		#pragma omp parallel
		{
			int id = omp_get_thread_num();
			if(id==0) { realthreads = omp_get_num_threads(); std::cout<<"\nRunning threads: "<<realthreads<<std::endl; }
			#pragma omp for schedule(dynamic)
				for(int p=0; p<Constants::population; p++)     
				{
					//std::cout<<"\nBouncing particle "<<p<<" "<<id<<std::flush;
					//Void Function for particle's motion. Involves RK4 for Nsteps. 
					//Detected particles are saved in ODPT object, which is passed here by reference.
					if( !(s1.compare(argv[1])) ) 		wpi(p, dstr[p], ODPT); 		//BELL + THE WAVE IS EVERYWHERE
					else if( !(s2.compare(argv[1])) )   wpi_ray(p, dstr[p], ODPT);  //LI   + RAY TRACING
				}
		}	
		std::cout<<"\n"<<"Joined"<<std::endl;
		real time2 = omp_get_wtime()-wtime ;
		std::cout<<"\nExecution time using "<<realthreads<<" thread(s), is: "<<time2<<std::endl;
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
		if(!dstr[p].id.empty()) //Only precipitated were saved, other vectors in the struct are empty.
		{
			precip_id.push_back(dstr[p].id.at(0));
			precip_lamda.push_back(dstr[p].lamda.at(0)); 
			precip_alpha.push_back(dstr[p].alpha.at(0));
			precip_aeq.push_back(dstr[p].aeq.at(0));
			precip_time.push_back(dstr[p].time.at(0));
		}
	}

	h5::File file("h5files/both_5000p_10s.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	//Detected particles
	h5::DataSet detected_lamda      = file.createDataSet("ODPT.lamda", ODPT.lamda);
	h5::DataSet detected_time       = file.createDataSet("ODPT.time", ODPT.time);
	h5::DataSet detected_id         = file.createDataSet("ODPT.id", ODPT.id);
	h5::DataSet detected_alpha      = file.createDataSet("ODPT.alpha", ODPT.alpha);

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


