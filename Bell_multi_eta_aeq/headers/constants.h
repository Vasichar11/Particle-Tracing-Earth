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
//e notations result in pico errors between py-c++.

//Define "constant" Parameters.

//--Particle initial parameters--//
const real L_shell=5; 			 				//L_shell of particle. Constant for now.
//const real aeq0_deg=20; //0-180 		     	//Initial equatorial pitch angle.
//const real aeq0=aeq0_deg*D2R;
const real lamda0_deg=10; //-30						
const real lamda0=lamda0_deg*D2R;	     		//Initial latitude, rads.
const real theta0_deg=0.001;            		//Initial wave normal angle.
const real theta0=theta0_deg*D2R;       

//---Wave parameters---//
const real By_wave=10*pow(10,-12);       		//10pT	         
const real Ekev0=500;               			//Initial energy keV
const real f_wave=2000; 			    		//Wave frequency in Hz. 2kHz
const real w_wave=2*M_PI*f_wave;        		//Wave angular frequency.

const real m_res=1;                           	//WPI resonance number (0=Landau resonance, 1= normal, counter-streaming resonance.)

//--Simulation time--//
const real t=1;        				  			//Simulation duration in seconds.
const real h=0.00001;							//Runge kutta stepsize.
const int64_t Nsteps=t/h; 			  			//Number of simulation steps, signed 8 byte integer.

const int64_t eta_dstr    = 3;					//How many sets of particles with different initial eta?											 
const int64_t aeq_dstr    = 3;					//Constituted by how many particles with different initial aeq?
const int64_t population = eta_dstr * aeq_dstr ;

//------------------------------ Eta distribution range and step. --------------------------------//
const real eta_start_d    = 0;     						 				//Degrees 	 
const real eta_end_d      = 360;   						 
const real eta_mid_d 	  = (eta_start_d + eta_end_d)/2;
const real eta_step_d 	  = (eta_end_d - eta_start_d)/(eta_dstr-1); 	//Step in degrees.

//------------------------------ Aeq distribution range and step. --------------------------------//					 
const real aeq_start_d    = 5;     						 				//Degrees
const real aeq_end_d      = 25;   			
const real aeq_mid_d 	  = (aeq_start_d + aeq_end_d)/2;
const real aeq_step_d	  = (aeq_end_d - aeq_start_d)/(aeq_dstr-1); 	//Step in degrees.
		                          
};	

