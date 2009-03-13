/************************************************************************
  			componentcap.cpp - 


**************************************************************************/

/*



   n       ||      p
---|--->---||------|------
       i   ||

Unp = Vn-Vp

i = Cd(Unp)/dt is used for the KCL


*/



#include "componentcapmna.h"
#include "algo.h"

componentCAPMNA::componentCAPMNA(dataCAP *d)
:componentCAP(d){
}
componentCAPMNA::~componentCAPMNA(){
}


void componentCAPMNA::stampBeforeInvertion(){
  stamp();
}
void componentCAPMNA::stamp(){
  ACE_CHECK_IERROR(mU,"componentCAPMNA::stamp mU null");
  ACE_CHECK_IERROR(mTenEq,"componentCAPMNA::stamp mTenEq null");
  
  bool dyn=false;
  //KCL equations
  algo::spls->KCL(mData.nodePos)->mCoefs[mU->mDynIndex]-=mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mU->mDynIndex]+=mData.value;
  
  //U=Vi-Vj
  mTenEq->mCoefs[mU->mIndex]-=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;
}
void componentCAPMNA::addUnknowns(){
  addTensionUnknown();
}
void componentCAPMNA::addEquations(){
  addTensionEquation();
}
