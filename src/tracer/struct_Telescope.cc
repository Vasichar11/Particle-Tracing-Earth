#include"struct_Telescope.h"

// Constructor - Initialize position of satellite
Telescope::Telescope(real lat, real L) : latitude_deg(lat), L_shell(L) {}

void Telescope::resize(size_t size){
	// Resize the vectors to the initial size
	latitude.resize(size);
	alpha.resize(size);
	aeq.resize(size);
	eta.resize(size);
	time.resize(size);
	id.resize(size);
}

// Returns true if particle crossed satellite
bool Telescope::crossing(real p1_latitude, real p2_latitude, real p_L_shell){
	bool crossed = false;

	if(p_L_shell == L_shell) {
		if ((p1_latitude - latitude_deg) * (p2_latitude - latitude_deg) < 0) {
		// if particle was below satellite and then above, product would be negative
			crossed = true;
		}
	}
	return crossed;
}

// Store detected particles in satellite's memory
void Telescope::store(int id, real latitude, real aeq, real alpha, real time){
	this->id.push_back(id);
	this->latitude.push_back(latitude);      	
	this->alpha.push_back(alpha);	
	this->aeq.push_back(aeq);			
	this->time.push_back(time);
}
