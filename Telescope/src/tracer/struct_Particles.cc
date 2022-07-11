#include "headers/struct_Particles.h"

void Particles::initialize(real eta0, real aeq0, real lamda0, real Ekin0, real zeta0, real time0)
{
	
	int k;
	const real Beq0 = Bmag_dipole(0);   		    //Beq isn't always Beq0?
	real Blam0      = Bmag_dipole(lamda0);
	real salpha0    = sin(aeq0)*sqrt(Blam0/Beq0);  //salpha = sin(aeq)*sqrt(Blam/Beq)
	if(aeq0*Constants::R2D>90)   k=1;			   //Both k=1 and k=0 are valid for every particle!!(?). This basically defines if its upward or downward.
	else					     k=0;			   //This way we distribute them: half upwards, half downwards.
	//srand (time(NULL)); //random seed using clock of computer
	//int k = rand() % 2; //rand number, 0 or 1.
	//If aeq<90 then k=0 then alpha0<90 <--> particle is northward, ppar>0
	//If aeq>90 then k=1 then alpha0>90 <--> particle is southward, ppar<0
	alpha0 = pow(-1,k)*asin(salpha0)+k*M_PI;       // sinx = a => x=(-1)^k * asin(a) + k*pi
	
	
	//Find momentum from energy.
	real Ejoule0=1.602176487E-16*Ekin0; //Kev to Joule
	real gama0=(Ejoule0/(Constants::m_e*pow(Constants::c,2))) + 1;
	real speed0=sqrt( 1 - (1/pow(gama0,2)) ) * Constants::c;
	//std::cout<<"\nBouncing period estimation: "<< (4*Constants::L_shell*Constants::Re/speed0)*(1.3 - 0.5*sin(aeq0)); //[Orlova1,Shprits2,2011]

	//Assign initial state.
	this->lamda0 = lamda0;
	this->aeq0   = aeq0;
	this->eta0   = eta0;
	this->alpha0 = alpha0;
	this->upar0  = speed0*cos(alpha0) ;
	this->uper0  = speed0*sin(alpha0) ;
	this->pper0  = gama0*Constants::m_e*uper0 ;
	this->ppar0  = gama0*Constants::m_e*upar0 ;
	this->zeta0  = zeta0 ;
	this->M_adiabatic0 = (pper0*pper0) / (2*Constants::m_e*Blam0) ;	
	this->Ekin0  = ((gama0-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15;	 //Joule back to Kev
	this->time0  = time0 ;
	//Particle is trapped, and won't escape until conditions are met.
	this->trapped  = true;	
	this->escaped  = false;
	this->negative = false; //May it develop negative P.A
	this->nan      = false; //May it develop NaN state
	this->high     = false; //May it develop NaN state


}


//Member function to save precipitating particle states. No need for vectors, just one saved value.
void Particles::escaping_state(int id, real new_lamda, real new_aeq, real new_alpha, real new_time)
{	
	id_lost     =	id;
	lamda_lost  =	new_lamda;
	aeq_lost  	=	new_aeq	;
	alpha_lost  =	new_alpha;
	time_lost   =	new_time;
	//ppar_lost   =	new_ppar;
	//pper_lost   =	new_pper;
	//zeta_lost   =	new_zeta;
	//uper_lost   =	new_uper;
	//upar_lost   =	new_upar;
	//eta_lost  	 =	new_eta	;
	//Ekin_lost   =	new_Ekin;
	//M_adiabatic_lost   =	new_M_adiabatic	;
}


//Member function to save particles that develop negative aeq.
void Particles::negative_state(int id)
{	
	id_neg     =	id;
}

//Member function to save precipitating that develop higher than 180 aeq.
void Particles::high_state(int id)
{	
	id_high     =	id;
}

//Member function to save precipitating that develop higher than na.
void Particles::nan_state(int id)
{	
	id_nan    =	id;
}

//Member function to save particle states(if needed).
//If need many values use vectors.

void Particles::save_state( int p,  real min_detadt, real mindetadt_time, real max_dEkin, real maxEkin_time, real  max_dPA, real maxdPA_time)
{
	saved_id = p;
	saved_max_dEkin = max_dEkin;
	saved_maxEkin_time = maxEkin_time;
	saved_max_dPA = max_dPA;
	saved_maxdPA_time = maxdPA_time;
	saved_min_detadt = min_detadt;
	saved_mindetadt_time = mindetadt_time;
	//this->id.push_back(id);
	//this->lamda.push_back(new_lamda);      				
	//this->ppar.push_back(new_ppar);		 
	//this->pper.push_back(new_pper);
	//this->alpha.push_back(new_alpha);	
	//this->aeq.push_back(new_aeq);
	//this->time.push_back(new_time);
	//this->zeta.push_back(new_zeta);
	//this->upar.push_back(new_upar);
	//this->uper.push_back(new_uper);	
	//this->deta_dt.push_back(new_deta_dt);
	//this->eta.push_back(new_eta);			
	//this->M_adiabatic.push_back(new_M_adiabatic); 
	//this->Ekin.push_back(new_Ekin); 
}	
