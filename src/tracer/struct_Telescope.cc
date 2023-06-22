#include"struct_Telescope.h"

//Constructor. Initialize position.
Telescope::Telescope(real latitude, real L_shell)
{
	L_shell   = L_shell;
	latitude  = latitude;
}
			
bool Telescope::crossing(real p1_latitude, real p2_latitude, real p_L_shell)
{
	bool crossed = false;

	if(p_L_shell == L_shell)
	{
		if( (latitude-p1_latitude)*(latitude-p2_latitude) < 0 ) 
		{ //e.g if particle was below satellite and then above, product would be negative.
			crossed = true;
		}
	}

	return crossed;
}

void Telescope::store(int id, real latitude, real aeq, real alpha, real time)
{
	this->id.push_back(id);
	this->latitude.push_back(latitude);      	
	this->alpha.push_back(alpha);	
	this->aeq.push_back(aeq);			
	this->time.push_back(time);
	//this->upar.push_back(upar);
	//this->uper.push_back(uper);	     
	//this->eta.push_back(eta);			
}
