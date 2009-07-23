/************************************************************************
  			equation.cpp 
**************************************************************************/
#include "equation.h"
#include <iostream>
using namespace std;

equation::equation(){
  mLine=-1;
  mCoefs=(ACE_DOUBLE *)0;
  mIsDyn = false;
  mIsLinear = true;
  mAvailable=true;
  mSize = 0;

}

void equation::allocMemory(int nb)
 {
   if(nb<1){
     ACE_INTERNAL_ERROR("equation alloc with nb<1");
     return;
   }
   mSize = nb;
    mCoefs = (ACE_DOUBLE*)calloc(mSize,sizeof(ACE_DOUBLE));
 }

 void equation::print(){
   for(int i =0; i <mSize;i++)
     cout <<"\t"<<mCoefs[i];
   cout <<"\n";
 }

equation::~equation(){
  if (mCoefs)
    free(mCoefs);
}


