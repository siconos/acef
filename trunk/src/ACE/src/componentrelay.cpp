/************************************************************************
  			componentrelay.cpp 
**************************************************************************/

#include "componentrelay.h"
#include "algo.h"

componentRELAY::componentRELAY(dataCOMP *d)
:component_LINEAR_NS(){
  ACE_CHECK_IERROR(d,"componentRELAY::componentRELAY : Diode data null");
  mData =(*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mNodeS=mData.nodeOut;
  mName = mData.name;
  mOffset = mData.voffset;
  mV2 = mData.vmax;
  mV1 = mData.vmin;
  
  mDimlambda=2;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;
  mType = ACE_TYPE_RELAY;
}
componentRELAY::~componentRELAY(){;}

void componentRELAY::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
  mVns=algo::spls->addinZns(ACE_TYPE_U,this);
  mIndiceStartZns= mVns->mIndexInVector;
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}

void componentRELAY::addEquations(){
  mEquation=algo::spls->addVdEquation(mName);
}
void componentRELAY::stamp(){
  int coef=1;
  if (mV2>mV1)
    coef=-1;
  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mNodeS)->mCoefs[i]=1;
  //because ie = ie'=0

  //VD laws
  //Zns = Vs-V0 ..
  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mNodeS);
  mEquation->mCoefs[i]+=+1;
  i= algo::spls->getIndexUnknown(ACE_TYPE_V,0);
  mEquation->mCoefs[i]-=1;
  mEquation->mCoefs[mVns->mIndex]-=1;

  
  //   |Z1|    |W+|
  //Y =|  |  L=|  |
  //   |W-|    |Z2|

  

  //  |Z1| |coef*(Uns-Vplus)  | | 0          | |coef*Uns| |0 | |-coef*Vplus |
  //Y=|  |=|                  |=|            |+|        |+|  |+|            |
  //  |W-| |-Ue+ W+ + offset  | | -(Vp - Vn) | |0       | |W+| |offset      |

  //| 0          |
  //|            |
  //| -(Vp - Vn) |
  if (mNodePos >0){
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodePos-1,-1);
  }
  if (mNodeNeg >0){
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeNeg-1,1);
  }
  //|coef*Uns|
  //|        |
  //|0       |
  algo::spls->mD1zns->setValue(mIndiceStartLambda,mIndiceStartZns,coef);
  //|0 |
  //|  |
  //|W+|
  algo::spls->mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda,1);
  //|-coef*Vplus |
  //|            |
  //|offset      |
  algo::spls->mD1s->setValue(mIndiceStartLambda,-coef*mV2);
  algo::spls->mD1s->setValue(mIndiceStartLambda+1,mOffset);


  

  //Zns = Vmloins -coef*Z2
  algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+1,-coef);
  algo::spls->mC1s->setValue(mIndiceStartZns,mV1);
}
void componentRELAY::print(){
  component_LINEAR_NS::print();
  printf("NodeS %d offset %f V1 %f V2 %f\n",mNodeS,mOffset,mV1,mV2);
  
}
