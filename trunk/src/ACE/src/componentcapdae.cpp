/************************************************************************
  			componentcap.cpp - 


**************************************************************************/

/*



   n       ||      p
---|--->---||------|------
       i   ||

Unp = Vn-Vp

Cd(Unp)/dt = i


*/



#include "componentcapdae.h"
#include "algo.h"

componentCAPDAE::componentCAPDAE(dataCAP *d)
:componentCAP(d){
}
componentCAPDAE::~componentCAPDAE(){
}


void componentCAPDAE::stampBeforeInvertion(){
  stamp();
}
void componentCAPDAE::stamp(){
  ACE_CHECK_IERROR(mU,"componentCAPDAE::stamp mU null");
  ACE_CHECK_IERROR(mI,"componentCAPDAE::stamp mI null");
  ACE_CHECK_IERROR(mTenEq,"componentCAPDAE::stamp mTenEq null");
  ACE_CHECK_IERROR(mDynEquation,"componentCAPDAE::stamp mDynEquation null");
  
  bool dyn=false;
  //KCL equations
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]-=1;
  algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;

  //U=Vi-Vj
  mTenEq->mCoefs[mU->mIndex]-=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;
  
  //Cu'=I
  mDynEquation->mCoefs[mU->mDynIndex]+=mData.value;
  mDynEquation->mCoefs[mI->mIndex]+=1;
  
}
