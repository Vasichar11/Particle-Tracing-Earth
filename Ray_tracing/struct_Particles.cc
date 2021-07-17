#include "headers/struct_Particles.h"

//Struct for particle's state throughout the iterations

//Constructor for initial state, called once when object is created.
Particles::Particles(real first_lamda, real first_zeta, real first_upar, real first_uper, real first_ppar, real first_pper, real first_alpha, real first_aeq, real first_aeqsu, real first_mu_ad_li, real first_time)
{ 
	this->lamda.push_back(first_lamda);      		//Make alternative constructor.			
	this->zeta.push_back(first_zeta);
	this->upar.push_back(first_upar);
	this->uper.push_back(first_uper);	     
	this->ppar.push_back(first_ppar);		 
	this->pper.push_back(first_pper);
	this->alpha.push_back(first_alpha);
	this->aeq.push_back(first_aeq);
	this->aeqsu.push_back(first_aeq);		
	this->mu_ad_li.push_back(first_mu_ad_li); //Add 
	this->time.push_back(first_time);
	//this_>L_shell of particle will be added
}   

//Member function to push back initial gyrophases.
void Particles::set_eta(real eta_distribution)
{
	eta.push_back(eta_distribution);
}

//Member function to push_back new state
void Particles::update_state(real new_lamda, real new_zeta, real new_uper, real new_upar, real new_ppar, real new_pper, real new_alpha, real new_aeq, real new_aeqsu, real new_eta, real new_mu_ad_li, real new_time)
{
	this->lamda.push_back(new_lamda);      				
	this->zeta.push_back(new_zeta);
	this->upar.push_back(new_upar);
	this->uper.push_back(new_uper);	     
	this->ppar.push_back(new_ppar);		 
	this->pper.push_back(new_pper);
	this->alpha.push_back(new_alpha);
	this->aeq.push_back(new_aeq);
	this->aeqsu.push_back(new_aeqsu);
	this->eta.push_back(new_eta);			
	this->mu_ad_li.push_back(new_mu_ad_li);			
	this->time.push_back(new_time);
}
	 	
