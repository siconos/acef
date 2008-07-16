/************************************************************************
  			componentind.cpp 
**************************************************************************/


/*

  n             p

--|-->--/\/\/\--|--
     i

Unp = Vn - Vp

Ldi/dt = Unp
     
*/


#include "componentind.h"
#include "algo.h"
componentIND::componentIND(dataIND * d){
  if(!d)
    ACE_ERROR("Inductor data null");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  if(ACE_IS_NULL(mData.value))
    ACE_ERROR("Inductor null");
  mType = ACE_TYPE_IND;

}

void componentIND::addUnknowns(){
  mI=algo::spls->addinx(ACE_TYPE_I,this);
}
void componentIND::addEquations(){
  mDynEquation=algo::spls->addIndEquation();
}
void componentIND::stamp(){
  //Li'=U
  int i;
  mDynEquation->mCoefs[mI->mDynIndex]+=mData.value;
  i=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos);
  mDynEquation->mCoefs[i]+=-1;
  i=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  mDynEquation->mCoefs[i]+=1;
  //KCL
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]+=-1;
  algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;
}
componentIND::~componentIND(){
  
}
