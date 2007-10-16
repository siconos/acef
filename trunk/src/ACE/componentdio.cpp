/************************************************************************
  			componentdio.cpp 
**************************************************************************/


#include "componentdio.h"
#include "algo.h"

componentDIO::componentDIO(dataDIO *d)
:componentNLINEAR(){
  if(!d)
    ACE_ERROR("Diode data null");
  mData =(*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
}

void componentDIO::addUnknowns(){
  mI=algo::sls.addinZns(ACE_TYPE_I,this);
}
void componentDIO::stamp(){
  int i=mI->mIndex;
  algo::sls.KCL(mData.nodeNeg)->mCoefs[i]+=1;
  algo::sls.KCL(mData.nodePos)->mCoefs[i]+=-1;  
}
componentDIO::~componentDIO(){
  
}
