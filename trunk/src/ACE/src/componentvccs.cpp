/************************************************************************
  			componentvccs.cpp 
**************************************************************************/

#include "componentvccs.h"
#include "algo.h"

componentVCCS::~componentVCCS(){;}
componentVCCS::componentVCCS(dataVCCS *d)
:componentLINEAR(){
  if(!d)
    ACE_ERROR("VCCS no data");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mType = ACE_TYPE_VCCS;
  if (ACE_IS_NULL(mData.value))
    ACE_WARNING("VCCS null");
}
void componentVCCS::addUnknowns(){
}
void componentVCCS::addEquations(){
}
void componentVCCS::stamp(){

  int drivNeg= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeDriverNeg);
  int drivPos= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeDriverPos);

  algo::spls->KCL(mData.nodeNeg)->mCoefs[drivNeg]-=mData.coef;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[drivPos]+=mData.coef;

  algo::spls->KCL(mData.nodePos)->mCoefs[drivNeg]+=mData.coef;
  algo::spls->KCL(mData.nodePos)->mCoefs[drivPos]-=mData.coef;
}
void componentVCCS::print(){
  componentLINEAR::print();
  printf(" driver neg %d, driver pos %d, coef %f",mData.nodeDriverNeg,mData.nodeDriverPos,mData.coef);
}



