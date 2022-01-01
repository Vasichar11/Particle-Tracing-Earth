#include "headers/interpolate.h"
#include <cassert>
#include <gnu/stubs-64.h>
#include <limits>
#include <iostream>


template <typename T,typename std::enable_if<std::is_arithmetic<T>::value>::type * = nullptr>
bool getLeftNeighbor(std::vector<T> data , T target,size_t& lNeighbor){

   for (size_t i=0; i<data.size()-1; i++){
      if (data.at(i)<=target && target <= data.at(i+1) ){
         lNeighbor=i;
         return true;
      }
   }
   //No neighbor found
   return false;
}

template <typename T,typename std::enable_if<std::is_arithmetic<T>::value>::type * = nullptr>
std::vector <T> interpolate(std::vector <T> time_new, std::vector <T> timef, std::vector <T> vector_csv)  {

#ifdef DEBUG
   //sanity check in debug mode
   bool sizeMatch = time_new.size()==timef.size() && timef.size() == vector_csv.size();
   assert(sizeMatch);
#endif

   size_t totalSize=time_new.size();
   std::vector <T> vector_int;vector_int.resize(totalSize);
   for (size_t i=0; i<totalSize; i++){
      size_t lNeighIndex;
      
      if ( !getLeftNeighbor(timef,time_new[i], lNeighIndex)){     
         vector_int[i]=std::numeric_limits<real>::quiet_NaN();
         std::cerr<<"NaN added to vector because no neighbors where found at  "<<__FILE__<<" "<<__LINE__<<std::endl;
         continue;
      }
      size_t rNeighIndex = lNeighIndex+1;
      T retval = vector_csv[lNeighIndex]+((time_new[i]-timef[lNeighIndex])*(vector_csv[rNeighIndex]-vector_csv[lNeighIndex])/(timef[rNeighIndex]-timef[lNeighIndex]));   
      vector_int[i]=retval;
   }
   return vector_int; 
}

