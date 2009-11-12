/************************************************************************
  			componentvcvs_SE1.cpp 
**************************************************************************/

#include "componentvcvs_SE1.h"
#include "algo.h"

componentVCVS_SE1::~componentVCVS_SE1(){;}
componentVCVS_SE1::componentVCVS_SE1(dataVCVS *d)
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
void componentVCVS_SE1::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
}
void componentVCVS_SE1::addEquations(){
  mEquation=algo::spls->addVdEquation();
}
void componentVCVS_SE1::stamp(){
  ACE_CHECK_IERROR(mI && mEquation,"componentVCVS_SE1::stamp no mI or no mEquation!!");
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]-=1;
  algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;

  //  int i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  //  mEquation->mCoefs[i]+=-1;
  algo::spls->mNodes[mData.nodeNeg]->stampV(-1,mEquation->mCoefs+algo::spls->mDimx,mEquation->mCoefs+2*algo::spls->mDimx);
  
  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos);
  //  mEquation->mCoefs[i]+=1;
  algo::spls->mNodes[mData.nodePos]->stampV(1,mEquation->mCoefs+algo::spls->mDimx,mEquation->mCoefs+2*algo::spls->mDimx);

  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeDriverNeg);
  //  mEquation->mCoefs[i]+=mData.coef;
  algo::spls->mNodes[mData.nodeDriverNeg]->stampV(mData.coef,mEquation->mCoefs+algo::spls->mDimx,mEquation->mCoefs+2*algo::spls->mDimx);

  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeDriverPos);
  //  mEquation->mCoefs[i]-=mData.coef;
  algo::spls->mNodes[mData.nodeDriverPos]->stampV(-mData.coef,mEquation->mCoefs+algo::spls->mDimx,mEquation->mCoefs+2*algo::spls->mDimx);

}
void componentVCVS_SE1::print(){
  componentLINEAR::print();
  printf(" driver neg %d, driver pos %d, coef %f",mData.nodeDriverNeg,mData.nodeDriverPos,mData.coef);
}
