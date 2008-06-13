/************************************************************************
  			componentvcvs.cpp 
**************************************************************************/

#include "componentvcvs.h"
#include "algo.h"

componentVCVS::~componentVCVS(){;}
componentVCVS::componentVCVS(dataVCVS *d)
:componentLINEAR(){
  if(!d)
    ACE_ERROR("VCVS no data");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mType = ACE_TYPE_VCVS;
  if (ACE_IS_NULL(mData.value))
    ACE_WARNING("VCVS null");
}
void componentVCVS::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
}
void componentVCVS::addEquations(){
  mEquation=algo::spls->addVdEquation();
}
void componentVCVS::stamp(){
  ACE_CHECK_IERROR(mI && mEquation,"componentVCVS::stamp no mI or no mEquation!!");
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]-=1;
  algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;

  int i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  mEquation->mCoefs[i]+=-1;
  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos);
  mEquation->mCoefs[i]+=1;

  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeDriverNeg);
  mEquation->mCoefs[i]+=mData.coef;
  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeDriverPos);
  mEquation->mCoefs[i]-=mData.coef;

}
void componentVCVS::print(){
  componentLINEAR::print();
  printf(" driver neg %d, driver pos %d, coef %f",mData.nodeDriverNeg,mData.nodeDriverPos,mData.coef);
}
