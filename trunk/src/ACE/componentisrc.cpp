/************************************************************************
  			componentisrc.cpp 
**************************************************************************/

#include "componentisrc.h"
#include "algo.h"

componentISRC::~componentISRC(){;}

componentISRC::componentISRC(dataISRC *d)
:componentLINEAR(){
  if(!d)
    ACE_ERROR("ISRC no data");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mCurrentValue=0;

  mType = ACE_TYPE_ISRC;
  if (ACE_IS_NULL(mData.value))
    ACE_WARNING("ISRC null");
}
void componentISRC::stamp(){
  //nothing, time dependant.
  //KCL
  //  algo::sls.KCL(mData.nodeNeg)->mCoefs[algo::sls.mRS]-=mData.value;
  //algo::sls.KCL(mData.nodePos)->mCoefs[algo::sls.mRS]+=mData.value;
}


void componentISRC::stampTimer(){
  ACE_DOUBLE newValue;
  ACE_CHECK_IERROR(getSourceValue("Isource",mData.id,&newValue),"componentISRC::stampTimer");
  //KCL
  algo::sls.KCL(mData.nodeNeg)->mCoefs[algo::sls.mRS]-=newValue - mCurrentValue;
  algo::sls.KCL(mData.nodePos)->mCoefs[algo::sls.mRS]+=newValue - mCurrentValue;

  mCurrentValue = newValue;
}



void componentISRC::print(){
  componentLINEAR::print();
  printf("\t value: %f\n",mData.value);
}
