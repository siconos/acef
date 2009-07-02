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
  //  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->mRS]-=mData.value;
  //algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->mRS]+=mData.value;
}


void componentISRC::stampTimer(){
  ACE_DOUBLE newValue;
  ACE_times[ACE_TIMER_TEST_7].start();
  // ParserGetSourceValue("Isource",mData.id,&newValue);
  getISRCValue(mData.id,&newValue);

  //KCL
  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->mRS]-=newValue - mCurrentValue;
  algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->mRS]+=newValue - mCurrentValue;
  mCurrentValue = newValue;
  ACE_times[ACE_TIMER_TEST_7].stop();
}



void componentISRC::print(){
  componentLINEAR::print();
  printf("\t value: %f\n",mData.value);
}
