/************************************************************************
  			componentdio.cpp 

Convention:
current go from Neg node to Pos node ==> current =<0
tension = V_neg - V_pos ==> tension >=0
stamp() to write : Ax'=Bx+CZs+DZns+s
                    0 =Ex+FZs+DZns+t
**************************************************************************/

#include "componentdio.h"
#include "algo.h"



componentDIO::componentDIO(dataDIO *d)
:componentNLINEAR(){
  ACE_CHECK_IERROR(d,"componentDIO::componentDIO : Diode data null");
  mData =(*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  ACE_CHECK_ERROR(mNodeNeg>=0 && mNodePos>=0 && mNodeNeg != mNodePos,"componentDIO::componentDIO");
  mDimlambda=1;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;
  mType = ACE_TYPE_DIO;
  

  
}

void componentDIO::addUnknowns(){
  mI=algo::spls->addinZns(ACE_TYPE_I,this);
  mIndiceStartZns= mI->mIndexInVector;
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}
void componentDIO::stamp(){
  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mData.nodeNeg)->mCoefs[i]-=1;
  algo::spls->KCL(mData.nodePos)->mCoefs[i]+=1;

  //stamp C1 and D1 system:
  algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,-1);
  if (mNodeNeg)
    algo::spls->mD1zs->setValue(mIndiceStartZns,mNodeNeg-1,1);
  if (mNodePos)
    algo::spls->mD1zs->setValue(mIndiceStartZns,mNodePos-1,-1);
  
}

componentDIO::~componentDIO(){
  
}
