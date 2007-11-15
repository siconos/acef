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

  mType = ACE_TYPE_ISRC;
  if (ACE_IS_NULL(mData.value))
    ACE_WARNING("ISRC null");
}
void componentISRC::stamp(){
  //KCL
  algo::sls.KCL(mData.nodeNeg)->mCoefs[algo::sls.mRS]-=mData.value;
  algo::sls.KCL(mData.nodePos)->mCoefs[algo::sls.mRS]+=mData.value;
}

void componentISRC::print(){
  componentLINEAR::print();
  printf("\t value: %f\n",mData.value);
}
