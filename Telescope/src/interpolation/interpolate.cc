#include "headers/interpolate.h"

std::vector <real> interpolate(std::vector <real> time_new, std::vector <real> timef, std::vector <real> vector_csv)  
{
    std::vector <real> vector_int;

    int size = time_new.size();
    vector_int.reserve(size);
    
    int x=1;
    real yp;
	for(size_t p=1; p<time_new.size(); p++)                        
    {             
		if(timef[x-1]<=time_new[p] && time_new[p]<=timef[x])      
        {   //Independent variable: time_new -> x
            //x1<=xp<=x2
            //x1,y1 coordinates Below the known xp value
            //x2,y2 coordinates Above the known xp value  
            //yp = y1 + (xp - x1)(y2-y1)/(x2-x1)                                                 

            //std::cout<< std::setprecision(16)<<std::scientific<<"\nx"<<x<<"\n csv 1 " << vector_csv[x-1) << "\n time_new p " << time_new[p) << "\n timef 1 " << timef[x-1) << "\n csv 2 " << vector_csv[x) <<"\n csv 1 "<< vector_csv[x-1)<<"\n timef 2 "<<timef[x)<<"\n timef 1 "<<timef[x-1) ; 
            
            yp = vector_csv[x-1] + ((time_new[p] - timef[x-1])*(vector_csv[x] - vector_csv[x-1])/(timef[x]-timef[x-1]));   
            vector_int.push_back(yp);
   		}
   		else
        {
            x++; //Change range of interpolation
            p--; //To loop again and find interpolant in this range. Has to be another way too...
        }
    }

    vector_int.insert(vector_int.begin(),vector_csv[0]); //Insert first value at start.

    return vector_int; 
}