/************************************************************************
  			componentarb.cpp 
**************************************************************************/

#include "componentarb.h"
#include "algo.h"

componentARB::~componentARB(){;}
componentARB::componentARB(dataARB *d)
:componentLINEAR(){
  if(!d)
    ACE_ERROR("ARB no data");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mType = ACE_TYPE_ARB;
  ACE_MESSAGE("B.... not managed");
}
void componentARB::addUnknowns(){
}
void componentARB::addEquations(){
}
void componentARB::stamp(){
  ACE_DOUBLE newValue;
}



