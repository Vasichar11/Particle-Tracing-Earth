#include "headers/distributions.h"

 std::array<real, Constants::aeq_dstr/2> half_norm(real mean, real std)
 {
    std::random_device seed;         //random seed
    std::mt19937 generator(seed());  //PRNG initialized with seed
    double number;
    std::array<real,Constants::aeq_dstr/2> half_dstr;

	//"Half-normal" dstr for da
    for(int i=0;i<Constants::aeq_dstr/2;i++)
    {
        std::normal_distribution<double> distribution(mean, std); //If inside the loop -> different values, why?
        number = distribution(generator); 
        if(number>=0) half_dstr.at(i) = number;  //take positives close to the mean
        else i--; //sample again
    }
    
    std::sort(half_dstr.begin(), half_dstr.end()); //Sort the da[] low->high
    
    return half_dstr;
 }

//"Symmetrical distribution" by mirroring half normal and shifting it.
 std::array<real, Constants::aeq_dstr> symmetrical_dstr(std::array<real, Constants::aeq_dstr/2> half_dstr, real shift)
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

