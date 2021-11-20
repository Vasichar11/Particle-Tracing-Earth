#include "headers/distributions.h"

 std::array<real, Constants::aeq_dstr/2 + 1> half_norm(real halfmean, real std)
 {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(halfmean, std); 
    double number;
    std::array<real,Constants::aeq_dstr/2 + 1> half_dstr;

	//"Half-normal" dstr for da
    for(int i=0;i<Constants::aeq_dstr/2;i++)
    {
        number = distribution(generator); 
        if(number>=0) half_dstr.at(i) = number; 
        else i--; //sample again
    }
    
    std::sort(half_dstr.begin(), half_dstr.end()); //Sort the da[] low->high
    
    return half_dstr;
 }



//"Symmetrical distribution" by mirroring half normal and shifting it.
 std::array<real, Constants::aeq_dstr> symmetrical_dstr(std::array<real, Constants::aeq_dstr/2 + 1> half_dstr, real shift)
{
    std::array<real, Constants::aeq_dstr> dstr;

	//Shift degrees. Symmetrical dstr arround mean=(shifting value)
	for(int i=0; i<Constants::aeq_dstr/2; i++)
    {
        dstr.at(Constants::aeq_dstr/2 + i)      = shift + half_dstr.at(i); //from half to end, i.e 50-99, 50 elements
        dstr.at(Constants::aeq_dstr/2 - i - 1)  = shift - half_dstr.at(i); //from half-1 to start, i.e 49-0, 50 elements
    }
    
    return dstr;
}

