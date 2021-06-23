typedef double real;  //Clarifies quantities
namespace Constants{

//Define Universal Costants.

real Re=6378137;  				        //Earth mean radius in m.
real c=2.997956376932163e8; 	        //Speed of light in m/s.
real m_e=9.10938291e-31;                //Electron mass in kg.
real q_e=1.602176565e-19;               //Electron charge in C.
real q_i=1.602176565e-19;               //Electron charge in C.
real m_O=2.67616E-026;                  //Oxygen mass in kg.    
real m_H=1.6726E-027;                   //Hydrogen mass in kg.  
real m_He=6.6904E-027;                  //Helium mass in kg.    
real e_el=5.105396765648739e5;	  
real eps0=8.854187817e-12;
real mu_0=M_PI*4*pow(10,-7);
real D2R=M_PI/180;
real R2D=1/D2R;
real B0=3.12*pow(10,-5);          		//Local uniform field in T.
real ne_0=10*pow(10,6);          		//m-3 //e notations result in pico errors between py-cpp.       

//Define "constant" Parameters.

//--Particle initial parameters--//
real L_shell=2; 				 		//L_shell of particle. Constant for now.
real aeq0_deg=70;  		     	 		//Initial equatorial pitch angle(defined as loss cone angle).
real aeq0=aeq0_deg*D2R;
real lamda0_deg=-2;//-9 				//resonance will occur at -5
real lamda0=lamda0_deg*D2R;	     		//Initial latitude, rads.
real theta0_deg=0.001;            		//Initial wave normal angle.
real theta0=theta0_deg*D2R;       

//---Wave parameters---//
real By_wave=1*pow(10,-12);       		//Byw=1pT.        
real Bw_tot_li=1.4142135*pow(10,-12);     
real Ekev0=511.3; 			            //Initial energy keV.
real f_wave=2500; 			   		    //Wave frequency in Hz.
real w_wave=2*M_PI*f_wave;        		//Wave angular frequency.

real m_res=1;                           //WPI resonance number (0=Landau resonance).
real pwr = pow(10,-4);                  //Poynting flux [W/m 2].

//--Simulation time--//
real t=1;        				  		//Simulation duration in seconds.
real h=0.00001;							//Runge kutta stepsize.
real pulse_duration=0.1;          		//Wave pulse duration in seconds.
real puls_dur=int(pulse_duration/h);	//Wave pulse duration in stepsize.
int64_t Nsteps=t/h; 			  		//Number of simulation steps, signed 8 byte integer.

//--Particle distribution initializations--//
double population = 6;		//6,24
double start_d    = 0;      //degrees 
double end_d      = 361;                                 
double lin_step_d = (end_d - start_d)/(population-1);  //Linearly spaced                                               

};	
