#pragma once
#include "common.h"
#include <cmath>

namespace Constants{

//--Universal parameters--//
const real Re=6378137;  				        //Earth mean radius in m.
const real c=2.997956376932163e8; 	        	//Speed of light in m/s.
const real m_e=9.10938291e-31;                	//Electron mass in kg.
const real q_e=1.602176565e-19;               	//Electron charge in C.
const real q_i=1.602176565e-19;               	//Electron charge in C.
const real m_O=2.67616E-026; 			     	//Oxygen mass in kg.    
const real m_H=1.6726E-027;                  	//Hydrogen mass in kg.  
const real m_He=6.6904E-027;                 	//Helium mass in kg.    
const real e_el=5.105396765648739e5;	  
const real eps0=8.854187817e-12;
const real mu_0=M_PI*4*pow(10,-7);
const real D2R=M_PI/180;
const real R2D=1/D2R;
const real B0=3.12*pow(10,-5);          		//Mean value of the field on the equator at the Earth’s surface.
const real ne_0=3*pow(10,6);                  	// 10/cm-3 => 10*10^6m/m-3



//--Constant parameters--//
//particle initials
const real L_shell=5; 			 				//L_shell of particle. Constant for now.
const real Ekin0=589.999;               			//Initial energy keV
const real eta0_deg=30;	   		  				//Initial particle phase (angle between Vperp and BwR).
const real eta0=eta0_deg*D2R;
const real aeq0_deg=69.844;                         //If single values are used.
const real aeq0=aeq0_deg*D2R;
const real lamda0_deg=-9;                       //If single values are used.
const real lamda0=lamda0_deg*D2R;
const real hm = 100*pow(10,3);                  //Minimum allowable mirroring altitude in m.
const real zm = (Re + hm)/(L_shell*Re);
const real alpha_lc = asin(sqrt(pow(zm,3)/sqrt(1+3*(1-zm)))); //Loss cone angle in radians.

//wave initials
const real f_wave=2000; 			    		//Wave frequency in Hz. 2kHz
const real w_wave=2*M_PI*f_wave;        		//Wave angular frequency.
const real m_res=1;                             //WPI resonance number (0=Landau resonance, 1= normal, counter-streaming resonance.)
const real theta0_deg=0.001;            		//Initial wave normal angle.
const real theta0=theta0_deg*D2R;   
//Li && ray tracing
const real pwr = pow(10,0);                     //From the power we get the intensity of the ray. Poynting flux [W/m 2].
const real pulse_duration=0.1;          		//Wave pulse duration in seconds.

//Bell code ONLY. Wave is everywhere0
const real By_wave=1*pow(10,-9);  

//--Simulation parameters--// 
const real h=0.00001;						    //Runge kutta stepsize. Has to be much less than the particle's gyroperiod?
const real puls_dur=int(pulse_duration/h);		//Wave pulse duration in stepsize.

//--Satellite parameters--//
const real telescope_lamda = 0; 

//--Distribution parameters--//
const int64_t aeq_dstr    = 1;  	    		    
const int64_t lamda_dstr  = 1; //maybe only one latitude per aeq	???
const int64_t eta_dstr    = 1; 	    		    
const int64_t Ekin_dstr   = 1;  	    		    

//P.A dstr
const real aeq_start_d    = 1;     //Degrees 		  			 				
const real aeq_end_d      = 179;  
//Latitude dstr
const real lamda_start_d  =-90;    		  			 				
const real lamda_end_d    = 90;  	
//Eta dstr
const real eta_start_d    = 1;      
const real eta_end_d      = 359;   //Issue if eta=0 or 360?				 
//Ekin dstr
const real Ekin_start      = 100;   //keV 	 
const real Ekin_end        = 1350;   		

//Normal dstr arguments
const real mean_aeq       = 90;  //Mean of the normal dstr
const real mean_lamda     = 0; 
const real mean_eta       = 180; 
const real mean_Ekin      = 600; 
const real stdev_aeq      = 20;  //Standard deviation of the normal dstr
const real stdev_eta      = 20;  //Standard deviation of the normal dstr
const real stdev_Ekin     = 20;  //Standard deviation of the normal dstr
const real max_stdev_lamda= 20;  //Standard deviation when lamda domain is the biggest. Is propotional to lamda range. 

//Latitude domain range precision
const real lamda_domain_step = 0.001;
				 
//Final population
const int64_t population  = aeq_dstr * (lamda_dstr * (eta_dstr * Ekin_dstr));   



};	

