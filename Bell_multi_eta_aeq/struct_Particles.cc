#include "headers/struct_Particles.h"

//Definition of struct.


//Struct for particle's state throughout the iterations
				
//Constructor.
Particles::Particles()
{

}

//Member function to push back initial gyrophases.
void Particles::set_eta(real initial_eta)
{
	this->eta.push_back(initial_eta);
}

//Member function to push back initial equatorial P.A.
void Particles::set_aeq(real initial_aeq)
{
	this->aeq.push_back(initial_aeq);
	this->aeq2.push_back(initial_aeq);
}

//Member function to calculate initial P.A, speed and momentum.
void Particles::calculations0(real Beq0, real initial_lamda,real initial_zeta, real initial_time, real initial_M_adiabatic, real initial_aeq, real initial_Ekev)
{

	//Find pitch angle at lamda0.
	real Blam0=Bmag_dipole(initial_lamda);
	real salpha0=sin(initial_aeq)*sqrt(Blam0/Beq0); //(2.20) Bortnik thesis. Is not true when M_adiabatic is "broken".
	
	//Check if valid
	//std::cout<<"\nsalpha0: "<<salpha0<< " when aeq0: " << initial_aeq*Constants::R2D;
	
	if(salpha0>1 || salpha0<-1) { throw 99; } //Temp exceptions. Are these ok?
	else if(salpha0 == 0) { throw 98; }
	
	real alpha0=asin(salpha0);					 	
	//Find momentum from energy.
	real Ejoule0=1.602176487E-16*initial_Ekev;
	real gama0=(Ejoule0/(Constants::m_e*pow(Constants::c,2))) + 1;
	real speed0=sqrt( 1 - (1/pow(gama0,2)) ) * Constants::c;
	real upar0=speed0*cos(alpha0);
	real uper0=speed0*sin(alpha0);
	real pper0=gama0*Constants::m_e*uper0;  //IF aeq=0->alpha0=0->pper=0 PROBLEM IN Bell_param, a2.(Exception 98)
	real ppar0=gama0*Constants::m_e*upar0;

	//Assign initial state.

	this->alpha.push_back(alpha0);
	this->alpha2.push_back(alpha0);
	this->upar.push_back(upar0);
	this->uper.push_back(uper0);
	this->pper.push_back(pper0);
	this->ppar.push_back(ppar0);
	this->lamda.push_back(initial_lamda);
	this->zeta.push_back(initial_zeta);
	this->time.push_back(initial_time);
	this->M_adiabatic.push_back(initial_M_adiabatic);	
}

//Member function to push_back new state
void Particles::update_state(real new_lamda, real new_zeta, real new_uper, real new_upar, real new_ppar, real new_pper, real new_alpha, real new_alpha2, real new_aeq, real new_aeq2, real new_eta, real new_M_adiabatic, real new_time)
{
	this->lamda.push_back(new_lamda);      				
	this->zeta.push_back(new_zeta);
	this->upar.push_back(new_upar);
	this->uper.push_back(new_uper);	     
	this->ppar.push_back(new_ppar);		 
	this->pper.push_back(new_pper);
	this->alpha.push_back(new_alpha);
	this->alpha2.push_back(new_alpha2);		
	this->aeq.push_back(new_aeq);
	this->aeq2.push_back(new_aeq2);
	this->eta.push_back(new_eta);			
	this->M_adiabatic.push_back(new_M_adiabatic); 
	this->time.push_back(new_time);
}	
