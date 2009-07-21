/************************************************************************
  			componentcap.cpp - 


**************************************************************************/

/*



   n       ||      p
---|--->---||------|------
       i   ||

i = Cd(Vn-Vp)/dt is used for the KCL


*/



#include "componentcapmna_v.h"
#include "algo.h"

componentCAPMNA_V::componentCAPMNA_V(dataCAP *d)
:componentCAP(d){
}
componentCAPMNA_V::~componentCAPMNA_V(){
}


void componentCAPMNA_V::stampBeforeInvertion(){
  stamp();
}
void componentCAPMNA_V::stamp(){
  
  //KCL equations
  //  algo::spls->KCL(mData.nodePos)->mCoefs[mU->mDynIndex]-=mData.value;
  algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]-=mData.value;
  algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodePos)]+=mData.value;
 
  //algo::spls->KCL(mData.nodeNeg)->mCoefs[mU->mDynIndex]+=mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=mData.value;

  //U=Vi-Vj
  /*mTenEq->mCoefs[mU->mIndex]-=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;*/
}
void componentCAPMNA_V::addUnknowns(){
  //  addTensionUnknown();
}
void componentCAPMNA_V::addEquations(){
  //  addTensionEquation();
}
