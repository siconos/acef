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

  
}

void componentDIO::addUnknowns(){
  mI=algo::sls.addinZns(ACE_TYPE_I,this);
  mIndiceStartZns= mI->mIndexInVector;
  mIndiceStartLambda= algo::sls.mDimLambda ;
  algo::sls.mDimLambda = algo::sls.mDimLambda + mDimlambda;
}
void componentDIO::stamp(){
  int i=mI->mIndex;
  //stamp equations.
  algo::sls.KCL(mData.nodeNeg)->mCoefs[i]-=1;
  algo::sls.KCL(mData.nodePos)->mCoefs[i]+=1;

  //stamp C1 and D1 system:
  algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,-1);
  if (mNodeNeg)
    algo::sls.mD1zs->setValue(mIndiceStartZns,mNodeNeg-1,1);
  if (mNodePos)
    algo::sls.mD1zs->setValue(mIndiceStartZns,mNodePos-1,-1);
  
}

componentDIO::~componentDIO(){
  
}
