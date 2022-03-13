#include "headers/struct_Particles.h"

void Particles::initialize(real eta0, real aeq0, real alpha0, real lamda0, real Ekev0, real Blam0, real zeta0, real time0, real lamda_start_d, real lamda_end_d)
{
	//Find momentum from energy.
	real Ejoule0=1.602176487E-16*Ekev0; //Kev to Joule
	real gama0=(Ejoule0/(Constants::m_e*pow(Constants::c,2))) + 1;
	real speed0=sqrt( 1 - (1/pow(gama0,2)) ) * Constants::c;
	//std::cout<<"\nBouncing period estimation: "<< (4*Constants::L_shell*Constants::Re/speed0)*(1.3 - 0.5*sin(aeq0)); //[Orlova1,Shprits2,2011]

	//Assign initial state.
	lamda_init = lamda0;
	aeq_init = aeq0;
	eta_init = eta0;
	alpha_init = alpha0;
	upar_init = speed0*cos(alpha0) ;
	uper_init = speed0*sin(alpha0) ;
	pper_init = gama0*Constants::m_e*uper_init ;
	ppar_init = gama0*Constants::m_e*upar_init ;
	zeta_init = zeta0 ;
	M_adiabatic_init = (pper_init*pper_init) / (2*Constants::m_e*Blam0) ;	
	Ekin_init = ((gama0-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15;	 //Joule back to Kev
	time_init = time0 ;
	escaped = false;	

	lamda_start_d = lamda_start_d;
	lamda_end_d   = lamda_end_d;


	//Else Invalid
	/*-1<salpha || salpha>1 => domain error
	aeq0=0||180 => salpha0 = 0 => alpha0 = 0 => pper0 = 0 => PROBLEM IN Bell_param, a2(division with 0).
	sin(pi) is actually not zero in C++... */

}


//Member function to save particle states
void Particles::save_state(int id, real new_lamda, real new_alpha, real new_aeq, real new_time)
{	
	id_lost     =	id;
	lamda_lost  =	new_lamda;
	alpha_lost  =	new_alpha;
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
	

	
	
	//If vectors are needed

	//this->id.push_back(id);
	//this->lamda.push_back(new_lamda);      				
	//this->zeta.push_back(new_zeta);
	//this->upar.push_back(new_upar);
	//this->uper.push_back(new_uper);	     
	//this->ppar.push_back(new_ppar);		 
	//this->pper.push_back(new_pper);
	//this->alpha.push_back(new_alpha);	
	//this->aeq.push_back(new_aeq);
	//this->deta_dt.push_back(new_deta_dt);
	//this->eta.push_back(new_eta);			
	//this->M_adiabatic.push_back(new_M_adiabatic); 
	//this->Ekin.push_back(new_Ekin); 
	//this->time.push_back(new_time);
}	
