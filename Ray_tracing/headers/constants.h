#pragma once
#include "common.h"
#include <cmath>

namespace Constants{

//Define Universal Costants.

const real Re=6378137;  				      //Earth mean radius in m.
const real c=2.997956376932163e8; 	          //Speed of light in m/s.
const real m_e=9.10938291e-31;                //Electron mass in kg.
const real q_e=1.602176565e-19;               //Electron charge in C.
const real q_i=1.602176565e-19;               //Electron charge in C.
const real m_O=2.67616E-026;                  //Oxygen mass in kg.    
const real m_H=1.6726E-027;                   //Hydrogen mass in kg.  
const real m_He=6.6904E-027;                  //Helium mass in kg.    
const real e_el=5.105396765648739e5;	  
const real eps0=8.854187817e-12;
const real mu_0=M_PI*4*pow(10,-7);
const real D2R=M_PI/180;
const real R2D=1/D2R;
const real B0=3.12*pow(10,-5);          		//Local uniform field in T.
const real ne_0=10*pow(10,6);          			//m-3 //e notations result in pico errors between py-cpp.       

//Define "constant" Parameters.

//--Particle initial parameters--//
const real L_shell=2; 				 			//L_shell of particle. Constant for now.
const real aeq0_deg=70;  		     	 		//Initial equatorial pitch angle(defined as loss cone angle).
const real aeq0=aeq0_deg*D2R;
const real lamda0_deg=-2;//-9 					//resonance will occur at -5
const real lamda0=lamda0_deg*D2R;	     		//Initial latitude, rads.
const real theta0_deg=0.001;            		//Initial wave normal angle.
const real theta0=theta0_deg*D2R;       

//---Wave parameters---//
const real By_wave=1*pow(10,-12);       		//Byw=1pT.        
const real Bw_tot_li=1.4142135*pow(10,-12);     
const real Ekev0=511.3; 			            //Initial energy keV.
const real f_wave=2500; 			   		    //Wave frequency in Hz.
const real w_wave=2*M_PI*f_wave;        		//Wave angular frequency.

const real m_res=1;                             //WPI resonance number (0=Landau resonance).
const real pwr = pow(10,-4);                    //Poynting flux [W/m 2].

//--Simulation time--//
const real t=1;        				  			//Simulation duration in seconds.
const real h=0.00001;							//Runge kutta stepsize.
const real pulse_duration=0.1;          		//Wave pulse duration in seconds.
const real puls_dur=int(pulse_duration/h);		//Wave pulse duration in stepsize.
const int Nsteps=t/h; 			  				//Number of simulation steps, signed 8 byte integer.

//--Particle distribution initializations--//
const real population = 6;		//6,24
const real start_d    = 0;      //degrees 
const real end_d      = 361;                                 
const real lin_step_d = (end_d - start_d)/(population-1);  //Linearly spaced                                               

};	
