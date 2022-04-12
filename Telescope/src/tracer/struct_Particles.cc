#include "headers/struct_Particles.h"

void Particles::initialize(real eta0, real aeq0, real alpha0, real lamda0, real Ekev0, real Blam0, real zeta0, real time0)
{
	//Find momentum from energy.
	real Ejoule0=1.602176487E-16*Ekev0; //Kev to Joule
	real gama0=(Ejoule0/(Constants::m_e*pow(Constants::c,2))) + 1;
	real speed0=sqrt( 1 - (1/pow(gama0,2)) ) * Constants::c;
	//std::cout<<"\nBouncing period estimation: "<< (4*Constants::L_shell*Constants::Re/speed0)*(1.3 - 0.5*sin(aeq0)); //[Orlova1,Shprits2,2011]
	
	//Assign initial state.
	this->lamda0 = lamda0;
	this->aeq0 = aeq0;
	this->eta0 = eta0;
	this->alpha0 = alpha0;
	this->upar0 = speed0*cos(alpha0) ;
	this->uper0 = speed0*sin(alpha0) ;
	this->pper0 = gama0*Constants::m_e*uper0 ;
	this->ppar0 = gama0*Constants::m_e*upar0 ;
	this->zeta0 = zeta0 ;
	this->M_adiabatic0 = (pper0*pper0) / (2*Constants::m_e*Blam0) ;	
	this->Ekin0 = ((gama0-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15;	 //Joule back to Kev
	this->time0 = time0 ;
	//Particle is trapped, and won't escape until conditions are met.
	this->trapped  = true;	
	this->escaped  = false;
	this->negative = false; //May it develop negative P.A
	this->nan      = false; //May it develop NaN state

	//Else Invalid
	/*-1<salpha || salpha>1 => domain error
	aeq0=0||180 => salpha0 = 0 => alpha0 = 0 => pper0 = 0 => PROBLEM IN Bell_param, a2(division with 0).
	sin(pi) is actually not zero in C++... */

}


//Member function to save precipitating particle states. No need for vectors, just one saved value.
void Particles::escaping_state(int id, real new_lamda, real new_aeq, real new_time)
{	
	id_lost     =	id;
	lamda_lost  =	new_lamda;
	aeq_lost  	 =	new_aeq	;
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


//Member function to save precipitating particle states. No need for vectors, just one saved value.
void Particles::negative_state(int id, real new_lamda, real new_aeq ,real new_alpha, real new_ppar, real new_pper, real new_time)
{	
	id_neg     =	id;
	lamda_neg  =	new_lamda;
	alpha_neg  =	new_alpha;
	aeq_neg    =	new_aeq	;
	time_neg   =	new_time;
	ppar_neg   =	new_ppar;
	pper_neg   =	new_pper;
	//zeta_neg   =	new_zeta;
	//uper_neg   =	new_uper;
	//upar_neg   =	new_upar;
	//eta_neg  	 =	new_eta	;
	//Ekin_neg   =	new_Ekin;
	//M_adiabatic_neg   =	new_M_adiabatic	;
}


//Member function to save all particle states(if needed). Need for vectors, save values in every step of the simulation.
/*
void Particles::save_state(int id, real new_lamda, real new_aeq, real new_ppar, real new_pper, real new_alpha, real new_time)
{
	this->id.push_back(id);
	this->lamda.push_back(new_lamda);      				
	this->ppar.push_back(new_ppar);		 
	this->pper.push_back(new_pper);
	this->alpha.push_back(new_alpha);	
	this->aeq.push_back(new_aeq);
	this->time.push_back(new_time);
	//this->zeta.push_back(new_zeta);
	//this->upar.push_back(new_upar);
	//this->uper.push_back(new_uper);	
	//this->deta_dt.push_back(new_deta_dt);
	//this->eta.push_back(new_eta);			
	//this->M_adiabatic.push_back(new_M_adiabatic); 
	//this->Ekin.push_back(new_Ekin); 
}	
*/