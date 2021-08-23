#include"headers/struct_Telescope.h"

//Constructor. Initialize position.
Telescope::Telescope(real lat, real L_parameter)
{
	L_shell   = L_parameter;
	latitude  = lat;
}
			

bool Telescope::crossing(real p1_lamda, real p2_lamda, real p_L_shell)
{
	bool crossed = false;

	if(p_L_shell == L_shell)
	{
		if( (latitude-p1_lamda)*(latitude-p2_lamda) < 0 ) 
		{ //e.g if particle was below satellite and then above, product would be negative.
			crossed = true;
		}
	}

	return crossed;
}

void Telescope::store(int id, real lamda, real zeta, real uper, real upar, real ppar, real pper, real alpha, real alpha2, real aeq, real aeq2, real eta, real M_adiabatic, real time)
{
	this->id.push_back(id);
	this->lamda.push_back(lamda);      				
	this->zeta.push_back(zeta);
	this->upar.push_back(upar);
	this->uper.push_back(uper);	     
	this->ppar.push_back(ppar);		 
	this->pper.push_back(pper);
	this->alpha.push_back(alpha);
	this->alpha2.push_back(alpha2);		
	this->aeq.push_back(aeq);
	this->aeq2.push_back(aeq2);
	this->eta.push_back(eta);			
	this->M_adiabatic.push_back(M_adiabatic); 
	this->time.push_back(time);
}


//void Telescope::update_look_dir()  		
//{	
//	angle  = sector * sector_range;					//Update look direction. 
//}		

//Post-processing
/*
for(sector=1; sector<=sector_num; sector++)
{
	look_dir();
	// (pa + 180) to find the supplementary angle, i.e the angle where the detector should be arround to count the particle.
	// B0 is reference.
	if( ((int)(pa + 180) % 360) >= (angle-sector_range/2)  &&  ((int)(pa + 180) % 360) <= (angle+sector_range/2) )	
	{	
		
	}	
} */
					