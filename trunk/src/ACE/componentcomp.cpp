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
  mVplus = mData.vmax;
  mVmoins = mData.vmin;
  ACE_CHECK_ERROR(mEpsilon >0,"componentCOMP::componentCOMP, epsilon==0");
  mD12= (mVmoins-mVplus)/(2.0*mEpsilon);
  mD11= -mD12;
  
  mDimlambda=2;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;
  mType = ACE_TYPE_COMP;
}
componentCOMP::~componentCOMP(){;}
void componentCOMP::addUnknowns(){
  mI=algo::sls.addinZs(ACE_TYPE_I,this);
  mVs=algo::sls.addinZns(ACE_TYPE_U,this);
  mIndiceStartZns= mVs->mIndexInVector;
  mIndiceStartLambda= algo::sls.mDimLambda ;
  algo::sls.mDimLambda = algo::sls.mDimLambda + mDimlambda;
}

void componentCOMP::addEquations(){
  mEquation=algo::sls.addVdEquation(mName);
}
void componentCOMP::stamp(){
  int i=mI->mIndex;
  //stamp equations.
  algo::sls.KCL(mNodeS)->mCoefs[i]=1;
  //because ie = ie'=0

  //Zns = Vj-Vk mB1..
  //VD laws
  i= algo::sls.getIndexUnknown(ACE_TYPE_V,mNodeS);
  mEquation->mCoefs[i]+=+1;
  i= algo::sls.getIndexUnknown(ACE_TYPE_V,0);
  mEquation->mCoefs[i]-=1;
  mEquation->mCoefs[mVs->mIndex]-=1;

  //Y=Vp-Vn + I*lambda +- epsilon
  //Vp-Vn
  if (mNodePos >0){
    algo::sls.mD1zs->setValue(mIndiceStartLambda,mNodePos-1,1);
    algo::sls.mD1zs->setValue(mIndiceStartLambda+1,mNodePos-1,1);
  }
  if (mNodeNeg >0){
    algo::sls.mD1zs->setValue(mIndiceStartLambda,mNodeNeg-1,-1);
    algo::sls.mD1zs->setValue(mIndiceStartLambda+1,mNodeNeg-1,-1);
  }
  //I*lambda
  algo::sls.mD1l->setValue(mIndiceStartLambda,mIndiceStartLambda,1);
  algo::sls.mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda+1,1);
  //+-epsilon
  algo::sls.mD1s->setValue(mIndiceStartLambda,mEpsilon);
  algo::sls.mD1s->setValue(mIndiceStartLambda+1,-mEpsilon);

  //Zns = Vplus +(d11 d12)lambda
  algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,mD11);
  algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+1,mD12);
  algo::sls.mC1s->setValue(mIndiceStartZns,mVplus);
}
void componentCOMP::print(){
  componentNLINEAR::print();
  printf("NodeS %d epsilon %f V+ %f V- %f\n",mNodeS,mEpsilon,mVplus,mVmoins);
  
}
