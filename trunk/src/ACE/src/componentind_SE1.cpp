/************************************************************************
  			componentind_SE1.cpp 
**************************************************************************/


/*

  n             p

--|-->--/\/\/\--|--
     i

Unp = Vn - Vp

Ldi/dt = Unp
     
*/


#include "componentind_SE1.h"
#include "algo.h"
componentIND_SE1::componentIND_SE1(dataIND * d){
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

void componentIND_SE1::addUnknowns(){
  mI=algo::spls->addinx(ACE_TYPE_I,this);
}
void componentIND_SE1::addEquations(){
  mDynEquation=algo::spls->addIndEquation();
}
void componentIND_SE1::stamp(){
  //Li'=U
  //int i;
  mDynEquation->mCoefs[mI->mDynIndex]+=mData.value;
  
  //i=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos);
  //  mDynEquation->mCoefs[i]+=-1;
  ACE_DOUBLE * pX = mDynEquation->mCoefs + algo::spls->mDimx;
  ACE_DOUBLE * pZs = mDynEquation->mCoefs + 2*algo::spls->mDimx;

  algo::spls->mNodes[mData.nodePos]->stampV(-1,pX,pZs);
  
  //i=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  //  mDynEquation->mCoefs[i]+=1;
  
  algo::spls->mNodes[mData.nodeNeg]->stampV(1,pX,pZs);

  //KCL
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]+=-1;
  algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;
}
componentIND_SE1::~componentIND_SE1(){
  
}
