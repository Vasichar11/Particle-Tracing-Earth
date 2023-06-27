#include"struct_Telescope.h"

// Constructor - Initialize position of satellite
Telescope::Telescope(real lat, real L) : latitude_deg(lat), L_shell(L) {}

// Particle position on the input - before and after position
// Returns true if particle crossed satellite
bool Telescope::crossing(real p1_latitude, real p2_latitude, real p_L_shell)
{
	bool crossed = false;

	if(p_L_shell == L_shell)
	{
		if ((p1_latitude - latitude_deg) * (p2_latitude - latitude_deg) < 0) 
		{ // if particle was below satellite and then above, product would be negative
			crossed = true;
		}
	}

	return crossed;
}

// Store detected particles in satellite's memory
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
