/************************************************************************
  			componentcap.cpp - 


**************************************************************************/


#include "componentcap.h"
#include "algo.h"

componentCAP::componentCAP(dataCAP *d)
:componentDYN(){
  if(!d)
    ACE_ERROR("Capacitor data null");    
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mType = ACE_TYPE_CAP;
  mU=0;
  mI=0;
  mDynEquation=0;
  if (ACE_IS_NULL(mData.value))
    ACE_ERROR("Capacitor null");
}
componentCAP::~componentCAP(){

}

void componentCAP::stampBeforeInvertion(){
  //KCL dyn equations
  if (mU){
    if (algo::sls.KCL(mData.nodeNeg)->mIsDyn){
      algo::sls.KCL(mData.nodeNeg)->mCoefs[mU->mDynIndex]=mData.value;
    }
    if (algo::sls.KCL(mData.nodePos)->mIsDyn){
      algo::sls.KCL(mData.nodePos)->mCoefs[mU->mDynIndex]=-mData.value;
    }
  }
  //Cu'=I
  if(mDynEquation){
    mDynEquation->mCoefs[mU->mDynIndex]=mData.value;
    mDynEquation->mCoefs[mI->mIndex]=1;
  }
}
void componentCAP::stampAfterInvertion(){
}

void componentCAP::addUnknowns(){
  if (!algo::sls.isUnknown(ACE_TYPE_U,this)){
       	mU=algo::sls.addinx(ACE_TYPE_U,this);
	
       	//graph.addEdge(n1,n2,mU);
  }
}
void componentCAP::addCurrentUnknown(){
  mI= algo::sls.addinZs(ACE_TYPE_I,this);
}
void componentCAP::addEquations(){
  mDynEquation=algo::sls.addCapEquation();
}
