#include "headers/is_in_packet.h"

//Check if the particle is within the wave's latitude.
//Return wave's index where there's matching latitude(WPI). Will be used for wave's data.
int is_in_packet(real min_lat, real max_lat, real lamda_tmp, int i, std::vector<real> &wave_lat)
{
	std::vector <real>::iterator ptr, requested_index;    //To find requested index
	int index;
	real diff, min_diff;                                  //To find minimum difference between latitudes

          
    if(min_lat < (lamda_tmp*Constants::R2D) && (lamda_tmp*Constants::R2D) < max_lat)
    {       
		min_diff = 100 ; 
		for(ptr = wave_lat.begin() + i; ptr < wave_lat.begin() + (Constants::puls_dur + i); ptr++)
        {
       	    diff = std::abs(*ptr - lamda_tmp*Constants::R2D) ;        //Return abs of difference
       	    if(diff < min_diff) 
       	    {          
			    requested_index = ptr; //Current requested index
			   	min_diff = diff;   	   //Current minimum                                                             
       	    }
       	}
	index = distance(wave_lat.begin(), requested_index); //Final requested index
	//std::cout<<"\nindex "<<int(index-Constants::puls_dur)<<"\n";
	}
	else
	{
		index = -1;	//Particle not in same latitude. No WPI.
	}

   
	return index ;
}
