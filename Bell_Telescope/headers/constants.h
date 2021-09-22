#pragma once
#include "common.h"
#include <cmath>
namespace Constants{

//Define Universal Costants.

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
const real ne_0=3*pow(10,6);                  	// 10/cm-3 => 10*10^6/m-3


//Define "constant" Parameters.

//--Particle initial parameters--//
const real L_shell=5; 			 				//L_shell of particle. Constant for now.
const real Ekev0=500;               			//Initial energy keV
const real aeq0_deg=40;            		     	//Initial equatorial pitch angle, deg.
const real aeq0=aeq0_deg*D2R;
const real lamda0_deg=0;						//Initial latitude, deg.
const real lamda0=lamda0_deg*D2R;	     		
const real eta0_deg=30;	   		  				//Initial particle phase (angle between Vperp and BwR).
const real eta0=eta0_deg*D2R;

//---Wave parameters---//
const real By_wave=10*pow(10,-12);       		//10pT	         
const real f_wave=2000; 			    		//Wave frequency in Hz. 2kHz
const real w_wave=2*M_PI*f_wave;        		//Wave angular frequency.
const real m_res=1;                           	//WPI resonance number (0=Landau resonance, 1= normal, counter-streaming resonance.)
const real theta0_deg=0.001;            		//Initial wave normal angle.
const real theta0=theta0_deg*D2R;   
const bool interaction = 0;						//0 for no interaction

//--Simulation parameters--//
const real t=5;        				  	    //Simulation duration in seconds.
const real h=0.00001;							//Runge kutta stepsize. Has to be much less than the particle's gyroperiod?
const int64_t Nsteps=t/h; 			  			//Number of simulation steps, signed 8 byte integer.

//--Satellite parameters--//
const real telescope_lamda = 30; //90 //0,10,20,30

//--Distribution parameters--//
const int64_t eta_dstr    = 1;					//Number of different values for each distribution.										 
const int64_t aeq_dstr    = 1000;  			    
const int64_t lamda_dstr  = 1;					
const int64_t population  = eta_dstr * aeq_dstr * lamda_dstr;

//------------------------------ Eta distribution range and step. --------------------------------//
const real eta_start_d    = eta0_deg;     //Degrees 	 
const real eta_end_d      = eta0_deg;   						 
const real eta_step_d     = 0;//(eta_end_d - eta_start_d)/(eta_dstr-1); 	//Step in degrees.

//------------------------------ Aeq distribution range and step. --------------------------------//					 
const real aeq_start_d    = 1;    					 				
const real aeq_end_d      = 179;  			
const real aeq_step_d	  = (aeq_end_d - aeq_start_d)/(aeq_dstr-1); 	
		                          
//----------------------------- Lamda distribution range and step. -------------------------------//					 
const real lamda_start_d  = lamda0_deg;     					 				
const real lamda_end_d    = lamda0_deg;   			
const real lamda_step_d	  = 0;//(lamda_end_d - lamda_start_d)/(lamda_dstr-1); 	
};	

