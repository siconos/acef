/************************************************************************
  			componentcap.cpp - 


**************************************************************************/


#include "componentcap.h"
#include "algo.h"

componentCAP::componentCAP(dataCAP *d)
:componentDYN(){
  ACE_CHECK_IERROR(d,"Capacitor data null");    
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name?mData.name:ACE_name;
  ACE_CHECK_WARNING(mData.name,"componentCAP::componentCAP : No name");
  mType = ACE_TYPE_CAP;
  mU=0;
  mI=0;
  mDynEquation=0;
  mTenEq=0;
  mICoefs=0;
  ACE_CHECK_IERROR(!ACE_IS_NULL(mData.value),"Capacitor null");
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
      mTenEq->mCoefs[mU->mIndex]-=1;
      mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
      mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;
    }
    if (algo::spls->KCL(mData.nodeNeg)->mIsDyn){
      if (mI)
	algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]-=1;
      else
	algo::spls->KCL(mData.nodeNeg)->mCoefs[mU->mDynIndex]+=mData.value;
      dyn = true;
    }else{
      ;//stamp KCL with current(afterinvertion)
    }
    if (algo::spls->KCL(mData.nodePos)->mIsDyn){
      if (mI)
	algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;
      else
	algo::spls->KCL(mData.nodePos)->mCoefs[mU->mDynIndex]-=mData.value;
      dyn = true;
    }else{
       ;//stamp KCL with current(afterinvertion)
   }
   ACE_CHECK_IWARNING(dyn,"componentCAP::stampBeforeInvertion : no dyn KCL");
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
  ACE_CHECK_IWARNING(mICoefs==0,"componentCAP::stampAfterInvertion alloc mIcoefs not null");
  mICoefs = (ACE_DOUBLE*)calloc(algo::spls->mNbUnknowns+1,sizeof(ACE_DOUBLE));
  algo::spls->getlinefromdxdt(mU->mDynIndex,mICoefs+algo::spls->mx.size());
  for (i=0;i<algo::spls->mNbUnknowns+1;i++){
      mICoefs[i]=mData.value*mICoefs[i];
  }

  printI();
  
  if (!algo::spls->KCL(mData.nodePos)->mIsDyn){
    for (i=0;i<algo::spls->mNbUnknowns+1;i++){
      algo::spls->KCL(mData.nodePos)->mCoefs[i]+=mICoefs[i];
    }
  }
  if (!algo::spls->KCL(mData.nodeNeg)->mIsDyn){
    for (i=0;i<algo::spls->mNbUnknowns+1;i++){
      algo::spls->KCL(mData.nodeNeg)->mCoefs[i]-=mICoefs[i];
    }
  }
}
void componentCAP::printI(){
  printf("print cap i coefs %s:\n",mName?mName:"no_name");
  if (mICoefs)
    for (int i =0; i <= algo::spls->mNbUnknowns; i++)
      printf("\t%f",mICoefs[i]);
  printf("\n");
}

void componentCAP::addTensionUnknown(){
  mU=algo::spls->addinx(ACE_TYPE_U,this);
}
void componentCAP::addCurrentUnknown(){
  mI= algo::spls->addinZs(ACE_TYPE_I,this);
}
void componentCAP::addCurrentEquation(){
  mDynEquation=algo::spls->addCapEquation();
}
void componentCAP::addTensionEquation(){
  mTenEq=algo::spls->addTenEquation();
}








void componentCAP::addUnknowns(){
    ACE_WARNING(" componentCAP::addUnknowns mustn't be called");
  }
void componentCAP::addEquations(){
     ACE_WARNING(" componentCAP::addEquations mustn't be called");
  }
