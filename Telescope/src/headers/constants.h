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
const real ne_0=3*pow(10,6);                  	// 10/cm-3 => 10*10^6/m-3


//--Constant parameters--//
//particle initials
const real L_shell=5; 			 				//L_shell of particle. Constant for now.
const real Ekev0=590;               			//Initial energy keV
const real aeq0_deg=170;            		    //Initial equatorial pitch angle, deg.
const real aeq0=aeq0_deg*D2R;
const real lamda0_deg=0;						//Initial latitude, deg.
const real lamda0=lamda0_deg*D2R;	     		
const real eta0_deg=30;	   		  				//Initial particle phase (angle between Vperp and BwR).
const real eta0=eta0_deg*D2R;

//wave initials
const real f_wave=2000; 			    		//Wave frequency in Hz. 2kHz
const real w_wave=2*M_PI*f_wave;        		//Wave angular frequency.
const real m_res=1;                           	//WPI resonance number (0=Landau resonance, 1= normal, counter-streaming resonance.)
const real theta0_deg=0.001;            		//Initial wave normal angle.
const real theta0=theta0_deg*D2R;   
//bell
const real By_wave=10*pow(10,-12);               
//li && ray tracing
const real pwr = pow(10,-4);                    //Poynting flux [W/m 2]. //for Ray tracing code  
const real pulse_duration=0.1;          		//Wave pulse duration in seconds.


//--Simulation parameters--//
const real t= 15;        				  	    //Simulation time in seconds.
const real t_nowpi = 10;                         //NoWPI time.
const real t_wpi = t - t_nowpi;                 //WPI time.
const real h=0.00001;						    //Runge kutta stepsize. Has to be much less than the particle's gyroperiod?
const int64_t Nsteps_wpi  = t_wpi/h; 			//WPI step count
const int64_t Nsteps_nowpi= t_nowpi/h;          //noWPI step count
const real puls_dur=int(pulse_duration/h);		//Wave pulse duration in stepsize.

//--Satellite parameters--//
const real telescope_lamda = 0; 

//--Distribution parameters--//
//const int64_t eta_dstr    = 1;					//Number of different values for each distribution.										 
const int64_t aeq_dstr    = 120;  			    
const int64_t lamda_dstr  = 100;					
const int64_t population  = lamda_dstr * aeq_dstr;   
//Aeq dstr					 
const real aeq_start_d    = 30;    		//Try narrower range because particles are breaking   			 				
const real aeq_end_d      = 150;  			
const real aeq_step_d	  = (aeq_end_d - aeq_start_d)/(aeq_dstr-1); 	//step for linspace(start,end,aeq_dstr) 
//Lamda dstr					 
const real lamda_start_d  = -90;			 				
const real lamda_end_d    = 90;	
//const real lamda_step_d	  =(lamda_end_d - lamda_start_d)/(lamda_dstr-1); 	

//Eta dstr
//const real eta_start_d    = eta0_deg;     //Degrees 	 
//const real eta_end_d      = eta0_deg;   						 
//const real eta_step_d     = 0;//(eta_end_d - eta_start_d)/(eta_dstr-1); 	//Step in degrees.
};	
