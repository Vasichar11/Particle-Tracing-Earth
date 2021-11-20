#include "headers/distributions.h"

//"Symmetrical distribution" by mirroring half normal and shifting it.
 std::array<real, Constants::aeq_dstr> symmetrical_dstr(real halfmean, real std, real shift)
{
    
    std::array<real, Constants::aeq_dstr> dstr;
    
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(halfmean, std); 
    double number;
    std::array<real,Constants::aeq_dstr/2 + 1> da;

	//"Half-normal" dstr for da
    for(int i=0;i<=Constants::aeq_dstr/2;i++)
    {
        number = distribution(generator); 
        if(number>=0) da.at(i) = number; 
        else i--; //sample again
    }
    std::sort(da.begin(), da.end()); //Sort the da[] low->high
  

    for(int i=0;i<Constants::aeq_dstr/2;i++)  std::cout<<"\n"<<da.at(i);

    std::cout<<"\n\n\n\n";
    
	//Shift degrees. Symmetrical dstr arround mean=(shifting value)
	for(int i=1; i<=Constants::aeq_dstr/2; i++)
    {
        dstr.at(Constants::aeq_dstr/2 + i - 1) = shift + da.at(i);
        dstr.at(Constants::aeq_dstr/2 - i)     = shift - da.at(i);
    }

    
    
    
    
    for(int i=0; i<Constants::aeq_dstr; i++)  std::cout<<"\n"<<dstr.at(i);
  
    return dstr;
}