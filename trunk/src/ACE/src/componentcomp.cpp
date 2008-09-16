/************************************************************************
  			componentcomp.cpp 
**************************************************************************/

#include "componentcomp.h"
#include "algo.h"

componentCOMP::componentCOMP(dataCOMP *d)
:componentNLINEAR(){
  ACE_CHECK_IERROR(d,"componentCOMP::componentCOMP : Diode data null");
  mData =(*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mNodeS=mData.nodeOut;
  mName = mData.name;
  mEpsilon = mData.vepsilon;
  mV2 = mData.vmax;
  mV1 = mData.vmin;
  ACE_CHECK_ERROR(mEpsilon >0,"componentCOMP::componentCOMP, epsilon==0");
  mD12= (mV1-mV2)/(2.0*mEpsilon);
  mD11= -mD12;
  
  mDimlambda=2;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;
  mType = ACE_TYPE_COMP;
}
componentCOMP::~componentCOMP(){;}
void componentCOMP::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
  mVns=algo::spls->addinZns(ACE_TYPE_U,this);
  mIndiceStartZns= mVns->mIndexInVector;
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}

void componentCOMP::addEquations(){
  mEquation=algo::spls->addVdEquation(mName);
}
void componentCOMP::stamp(){
  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mNodeS)->mCoefs[i]=1;
  //because ie = ie'=0

  //Zns = Vj-Vk mB1..
  //VD laws
  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mNodeS);
  mEquation->mCoefs[i]+=+1;
  i= algo::spls->getIndexUnknown(ACE_TYPE_V,0);
  mEquation->mCoefs[i]-=1;
  mEquation->mCoefs[mVns->mIndex]-=1;

  //Y=Vp-Vn + I*lambda +- epsilon
  //Vp-Vn
  if (mNodePos >0){
    algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodePos-1,1);
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodePos-1,1);
  }
  if (mNodeNeg >0){
    algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodeNeg-1,-1);
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeNeg-1,-1);
  }
  //I*lambda
  algo::spls->mD1l->setValue(mIndiceStartLambda,mIndiceStartLambda,1);
  algo::spls->mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda+1,1);
  //+-epsilon
  algo::spls->mD1s->setValue(mIndiceStartLambda,mEpsilon);
  algo::spls->mD1s->setValue(mIndiceStartLambda+1,-mEpsilon);

  //Zns = Vplus +(d11 d12)lambda
  algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,mD11);
  algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+1,mD12);
  algo::spls->mC1s->setValue(mIndiceStartZns,mV2);
}
void componentCOMP::print(){
  componentNLINEAR::print();
  printf("NodeS %d epsilon %f V1 %f V2 %f\n",mNodeS,mEpsilon,mV1,mV2);
  
}
