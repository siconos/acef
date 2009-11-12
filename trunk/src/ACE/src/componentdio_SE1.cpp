/************************************************************************
  			componentdio_SE1.cpp 

Convention:
current go from Neg node to Pos node ==> current =<0
tension = V_neg - V_pos ==> tension >=0
stamp() to write : Ax'=Bx+CZs+DZns+s
                    0 =Ex+FZs+DZns+t
**************************************************************************/

#include "componentdio_SE1.h"
#include "algo.h"



componentDIO_SE1::componentDIO_SE1(dataDIO *d)
:component_LINEAR_NS(){
  ACE_CHECK_IERROR(d,"componentDIO_SE1::componentDIO_SE1 : Diode data null");
  mData =(*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  ACE_CHECK_ERROR(mNodeNeg>=0 && mNodePos>=0 && mNodeNeg != mNodePos,"componentDIO_SE1::componentDIO_SE1");
  mDimlambda=1;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;
  mType = ACE_TYPE_DIO;
  mThreshold=ACE_DIODE_THRESHOLD;
  cout<<"*********************************************************************"<<endl;
  cout<<"*************************WARNING, DIODE THRESHOLD IS : "<<mThreshold<<endl;
  

  
}

void componentDIO_SE1::addUnknowns(){
  mI=algo::spls->addinZns(ACE_TYPE_I,this);
  mIndiceStartZns= mI->mIndexInVector;
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}
void componentDIO_SE1::stamp(){
  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mData.nodeNeg)->mCoefs[i]-=1;
  algo::spls->KCL(mData.nodePos)->mCoefs[i]+=1;

  //stamp C1 and D1 system:
  algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,-1);
  ACE_CHECK_IERROR(mIndiceStartLambda == mIndiceStartZns,"dio model component, mIndiceStartLambda != mIndiceStartZns");
  //Bug since the begining, previous code : A lot of test has been validated is the following wrong code :
  //   if (mNodeNeg)
  //     algo::spls->mD1zs->setValue(mIndiceStartZns,mNodeNeg-1,1);
  //   if (mNodePos)
  //     algo::spls->mD1zs->setValue(mIndiceStartZns,mNodePos-1,-1);

  //may be mIndiceStartLambda == mIndiceStartZns. Else big problem.
  if (ACE_FORMULATION !=  ACE_FORMULATION_MNA_V && ACE_FORMULATION !=  ACE_FORMULATION_STAMP_ONLY){
    if (mNodeNeg){
      //      algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodeNeg-1,1);
      algo::spls->mNodes[mNodeNeg]->stampV(1,mIndiceStartLambda,algo::spls->mD1x,algo::spls->mD1zs);
    }
    if (mNodePos){
      //      algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodePos-1,-1);
      algo::spls->mNodes[mNodePos]->stampV(-1,mIndiceStartLambda,algo::spls->mD1x,algo::spls->mD1zs);
    }
  }else{
    if (mNodeNeg)
      algo::spls->mD1x->setValue(mIndiceStartLambda,mNodeNeg-1,1);
    if (mNodePos)
      algo::spls->mD1x->setValue(mIndiceStartLambda,mNodePos-1,-1);
  }
  algo::spls->mD1s->setValue(mIndiceStartLambda,mThreshold);
  
}

componentDIO_SE1::~componentDIO_SE1(){
  
}
void componentDIO_SE1::print (){
  component_LINEAR_NS::print();
   printf("diode threshold : %lf \n",mThreshold);
  ;
}
