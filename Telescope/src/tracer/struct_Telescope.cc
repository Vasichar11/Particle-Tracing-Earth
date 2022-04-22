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

void Telescope::store(int id, real lamda, real aeq, real alpha, real time)
{
	this->id.push_back(id);
	this->lamda.push_back(lamda);      	
	this->alpha.push_back(alpha);	
	this->aeq.push_back(aeq);			
	this->time.push_back(time);
	//this->upar.push_back(upar);
	//this->uper.push_back(uper);	     
	//this->eta.push_back(eta);			
}
