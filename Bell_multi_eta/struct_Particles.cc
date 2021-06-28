#include "headers/struct_Particles.h"

//Definition of struct.


//Struct for particle's state throughout the iterations
				
//Constructor for initial state, called once when object is created.
Particles::Particles(real first_lamda, real first_zeta, real first_upar, real first_uper, real first_ppar, real first_pper, real first_alpha, real first_alpha2, real first_aeq, real first_aeq2,real first_M_adiabatic,real first_time)
{ 
	this->lamda.push_back(first_lamda);      			//Make alternative constructor.			
	this->zeta.push_back(first_zeta);
	this->upar.push_back(first_upar);
	this->uper.push_back(first_uper);	     
	this->ppar.push_back(first_ppar);		 
	this->pper.push_back(first_pper);
	this->alpha.push_back(first_alpha);
	this->alpha2.push_back(first_alpha2);
	this->aeq.push_back(first_aeq);
	this->aeq2.push_back(first_aeq);
	this->time.push_back(first_time);
	this->M_adiabatic.push_back(first_M_adiabatic);	//Add and M adiabatic in struct.
	//this->eta.push_back(first_eta);		
	//this->L_shell of particle will be added
}   
//Member function to push back initial gyrophases.
void Particles::set_eta(real eta_distribution)
{
	this->eta.push_back(eta_distribution);
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
