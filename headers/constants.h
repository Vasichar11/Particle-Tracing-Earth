#pragma once
#include "common.h"
#include <cmath>

//--Universal parameters--//
namespace Universal
{
    const real Re = 6378137; // Earth mean radius in m
    const real c = 2.997956376932163e8; // Speed of light in m/s
    const real m_e = 9.10938291e-31; // Electron mass in kg
    const real q_e = 1.602176565e-19; // Electron charge in C
    const real q_i = 1.602176565e-19; // Electron charge in C
    const real m_O = 2.67616E-026; // Oxygen mass in kg  
    const real m_H = 1.6726E-027; // Hydrogen mass in kg
    const real m_He = 6.6904E-027; // Helium mass in kg    
    const real e_el = 5.105396765648739e5;	  
    const real eps0 = 8.854187817e-12;
    const real mu_0 = M_PI*4*pow(10,-7);
    const real D2R = M_PI/180;
    const real R2D = 1/D2R;
    const real B0 = 3.12*pow(10,-5); // Mean value of the field on the equator at the Earth’s surface
    const real ne_0 = 3*pow(10,6); // 10/cm-3 => 10*10^6m/m-3
}

//--Wave parameters--// 
namespace Wave
{
    const real f_wave = 2000; // Wave frequency in Hz
    const real w_wave = 2*M_PI*f_wave; // Wave angular frequency
    const real m_res = 1; // WPI resonance number (0=Landau resonance, 1= normal, counter-streaming resonance.)

    //Only Li code. Ray tracing
    const real pwr = pow(10,0); // From the power we get the intensity of the ray. Poynting flux [W/m 2]
    const real pulse_duration = 0.1; // Wave pulse duration in seconds

    //Only Bell code. Wave is everywhere
    const real theta0_deg=0.001; // Initial wave norml angle
    const real theta0 = theta0_deg*Universal::D2R;   
    const real By_wave = 1*pow(10,-9);  
}

//--Distribution parameters--//
namespace Distribution
{
    const real L_shell = 2; // L_shell of particle. Constant for now
    const int64_t population = 100;   
}

// P.A dstr
namespace Aeq_dstr
{
    const real start_deg = 1;  		  			 				
    const real end_deg = 179;  
    const real mean = 90; // Mean of the normal dstr
    const real stdev = 20; // Standard deviation of the normal dstr
    const real steps = 10; // To distribute evenly
    const real value_deg = 37.6708693574307; // If single values are used
    const real value = value_deg * Universal::D2R;
}
    
// Latitude dstr
namespace Lat_dstr
{
    const real start_deg = -90;    		  			 				
    const real end_deg = 90;  	
    const real mean = 0; // Mean of the normal dstr
    const real max_stdev = 20; // Standard deviation when latitude domain is the biggest. Is propotional to latitude range 
    const real domain_step = 0.001; // Latitude domain range precision, degrees
    const real steps = 10; // To distribute evenly
    const real value_deg = 16.2782913046611; // If single values are used
    const real value = value_deg * Universal::D2R;
}
    
// Eta dstr
namespace Eta_dstr
{
    const real start_deg = 1;      
    const real end_deg = 359; // Issue if eta=0 or 360?				 
    const real mean = 180; // Mean of the normal dstr
    const real stdev = 20; // Standard deviation of the normal dstr
    const real steps = 10; // To distribute evenly
    const real value_deg = 5.01547617274519; // Initial particle phase (angle between Vperp and BwR)
    const real value = value_deg * Universal::D2R;
}
    
// Ekin dstr
namespace Ekin_dstr
{
    const real start = 100; // keV 	 
    const real end = 1350; // keV 
    const real mean = 600; // Mean of the normal dstr		
    const real stdev = 20; // Standard deviation of the normal dstr
    const real steps = 10; // To distribute evenly
    const real value = 850.508389025729; // Initial energy keV
}

//--Simulation parameters--// 
namespace Simulation
{
    const real hm = 100*pow(10,3); // Minimum allowable mirroring altitude in m
    const real zm = (Universal::Re + hm)/(Distribution::L_shell*Universal::Re);
    const real alpha_lc = asin(sqrt(pow(zm,3)/sqrt(1+3*(1-zm)))); // Loss cone angle in radians
    const real h = 0.00001; // Runge kutta stepsize. Should be much smaller than the particle's gyroperiod
    const real puls_dur = int(Wave::pulse_duration/h); // Wave pulse duration in stepsize
}

//--Satellite parameters--//
namespace Satellite
{
    const real latitude_deg = 0; // Latitude of the satellite
    const real L_shell = Distribution::L_shell; // L_shell of the satellite
}