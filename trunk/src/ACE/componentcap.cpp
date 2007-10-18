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
  mTenEq=0;
  mICoefs=0;
  if (ACE_IS_NULL(mData.value))
    ACE_ERROR("Capacitor null");
}
componentCAP::~componentCAP(){
  if (mICoefs)
    free (mICoefs);

}

void componentCAP::stampBeforeInvertion(){
  //KCL dyn equations
  if (mU){
    bool dyn=false;
    if (mTenEq){
      mTenEq->mCoefs[mU->mIndex]+=1;
      mTenEq->mCoefs[algo::sls.getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
      mTenEq->mCoefs[algo::sls.getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;
    }
    if (algo::sls.KCL(mData.nodeNeg)->mIsDyn){
      algo::sls.KCL(mData.nodeNeg)->mCoefs[mU->mDynIndex]+=mData.value;
      dyn = true;
    }else{
      ;//stamp KCL with current(afterinvertion)
    }
    if (algo::sls.KCL(mData.nodePos)->mIsDyn){
      algo::sls.KCL(mData.nodePos)->mCoefs[mU->mDynIndex]-=mData.value;
      dyn = true;
    }else{
       ;//stamp KCL with current(afterinvertion)
   }
    if (!dyn)
      ACE_INTERNAL_WARNING("componentCAP::stampBeforeInvertion : no dyn KCL");
  }else{
    ACE_INTERNAL_WARNING("componentCAP::stampBeforeInvertion component without U.");
  }
  //Cu'=I
  if(mDynEquation){
    mDynEquation->mCoefs[mU->mDynIndex]+=mData.value;
    mDynEquation->mCoefs[mI->mIndex]+=1;
  }
}
void componentCAP::stampAfterInvertion(){
  int i;
  if (!mICoefs){
    ACE_INTERNAL_WARNING("componentCAP::stampAfterInvertion mICoefs null.");
    return;
  }
  if (!algo::sls.KCL(mData.nodePos)->mIsDyn){
    for (i=0;i<algo::sls.mNbUnknowns;i++){
      algo::sls.KCL(mData.nodePos)->mCoefs[i]+=mICoefs[i];
    }
  }
  if (!algo::sls.KCL(mData.nodeNeg)->mIsDyn){
    for (i=0;i<algo::sls.mNbUnknowns;i++){
      algo::sls.KCL(mData.nodeNeg)->mCoefs[i]+=mICoefs[i];
    }
  }
}

void componentCAP::addTensionUnknown(){
  mU=algo::sls.addinx(ACE_TYPE_U,this);
}
void componentCAP::addCurrentUnknown(){
  mI= algo::sls.addinZs(ACE_TYPE_I,this);
}
void componentCAP::addCurrentEquation(){
  mDynEquation=algo::sls.addCapEquation();
}
void componentCAP::addTensionEquation(){
  mTenEq=algo::sls.addTenEquation();
}








void componentCAP::addUnknowns(){
    ACE_WARNING(" componentCAP::addUnknowns mustn't be called");
  }
void componentCAP::addEquations(){
     ACE_WARNING(" componentCAP::addEquations mustn't be called");
  }
