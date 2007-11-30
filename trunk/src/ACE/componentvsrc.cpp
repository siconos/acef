/************************************************************************
  			componentvsrc.cpp 
**************************************************************************/

#include "componentvsrc.h"
#include "algo.h"

componentVSRC::~componentVSRC(){;}
componentVSRC::componentVSRC(dataVSRC *d)
:componentLINEAR(){
  if(!d)
    ACE_ERROR("VSRC no data");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mType = ACE_TYPE_VSRC;
  mCurrentValue=0;
  if (ACE_IS_NULL(mData.value))
    ACE_WARNING("VSRC null");
}
void componentVSRC::addUnknowns(){
  mI=algo::sls.addinZs(ACE_TYPE_I,this);
}
void componentVSRC::addEquations(){
  mEquation=algo::sls.addVdEquation();
}
void componentVSRC::stamp(){
  ACE_DOUBLE newValue;

  //KCL
  int i;
  ACE_CHECK_IERROR(mI && mEquation,"componentVSRC::stamp no mI or no mEquation!!");
  algo::sls.KCL(mData.nodeNeg)->mCoefs[mI->mIndex]-=1;
  algo::sls.KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;
  
  //VD laws
  i= algo::sls.getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  mEquation->mCoefs[i]+=-1;
  i= algo::sls.getIndexUnknown(ACE_TYPE_V,mData.nodePos);
  mEquation->mCoefs[i]+=1;
}

void componentVSRC::stampTimer(){
  ACE_DOUBLE newValue;
  getSourceValue("Vsource",mData.id,&newValue);
  mEquation->mCoefs[algo::sls.mRS]-=newValue - mCurrentValue;

  mCurrentValue = newValue;
}


