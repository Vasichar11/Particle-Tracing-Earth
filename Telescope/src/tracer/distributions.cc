#include "headers/distributions.h"

//"Symmetrical distribution" by mirroring half normal and shifting it.
 std::array<real, Constants::aeq_dstr> symmetrical_dstr(real halfmean, real std, real shift)
{
    
    std::array<real, Constants::aeq_dstr> dstr;
    
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(halfmean, std); 
    double number;
    std::array<real,Constants::aeq_dstr> da;

	//"Half-normal" dstr for da
    for(int i=0;i<=Constants::aeq_dstr/2;i++)
    {
        number = distribution(generator); 
        if(number>=0) da.at(i) = number; 
        else i--; //sample again
    }
    std::cout<<"\n"<<da.begin();
    std::sort(da.begin(), da.begin()+Constants::aeq_dstr/2,std::greater<int>()); //Sort the da[] high->low
  
    for(int i=1; i<Constants::aeq_dstr/2; i++)  da.at(Constants::aeq_dstr/2 + i) =  da.at(Constants::aeq_dstr/2 - i);

    for(int i=0;i<Constants::aeq_dstr;i++)  std::cout<<"\n"<<da.at(i);


	//Shift degrees. Symmetrical dstr arround mean=(shifting value)
	for(int i=0; i<Constants::aeq_dstr; i++)  dstr.at(i) += shift;
  
    return dstr;
}