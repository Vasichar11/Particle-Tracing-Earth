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
void Particles::save_state(real new_lamda, real new_alpha, real new_aeq, real new_ppar, real new_pper, real new_time)
{	//define size
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
