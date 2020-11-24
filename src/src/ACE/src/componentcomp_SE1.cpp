/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file componentcomp_SE1.cpp

*/
/************************************************************************
  			componentcomp_SE1.cpp 
**************************************************************************/

#include "componentcomp_SE1.h"
#include "algo.h"

componentCOMP_SE1::componentCOMP_SE1(dataCOMP *d)
:component_LINEAR_NS(){
  ACE_CHECK_IERROR(d,"componentCOMP_SE1::componentCOMP_SE1 : Diode data null");
  mData =(*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mNodeS=mData.nodeOut;
  mName = mData.name;
  mEpsilon = mData.vepsilon;
  mV2 = mData.vmax;
  mV1 = mData.vmin;
  ACE_CHECK_ERROR(mEpsilon >0,"componentCOMP_SE1::componentCOMP_SE1, epsilon==0");
  mD12= (mV1-mV2)/(2.0*mEpsilon);
  mD11= -mD12;
  
  mDimlambda=2;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;
  mType = ACE_TYPE_COMP;
}
componentCOMP_SE1::~componentCOMP_SE1(){;}
void componentCOMP_SE1::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
  mVns=algo::spls->addinZns(ACE_TYPE_U,this);
  mIndiceStartZns= mVns->mIndexInVector;
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}

void componentCOMP_SE1::addEquations(){
  mEquation=algo::spls->addVdEquation(mName);
}
void componentCOMP_SE1::stamp(){

  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mNodeS)->mCoefs[i]=1;
  //because ie = ie'=0

  //Zns = Vj-Vk mB1..
  //VD laws
  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mNodeS);
  //  mEquation->mCoefs[i]+=+1;
  algo::spls->mNodes[mNodeS]->stampV(1,mEquation->mCoefs +algo::spls->mDimx,mEquation->mCoefs +2*algo::spls->mDimx);
    
  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,0);
  algo::spls->mNodes[0]->stampV(-1,mEquation->mCoefs +algo::spls->mDimx,mEquation->mCoefs +2*algo::spls->mDimx);
  //  mEquation->mCoefs[i]-=1;
  
  mEquation->mCoefs[mVns->mIndex]-=1;

  //Y=Vp-Vn + I*lambda +- epsilon
  //Vp-Vn
  if (mNodePos >0){
    //    algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodePos-1,1);
    algo::spls->mNodes[mNodePos]->stampV(1,mIndiceStartLambda,algo::spls->mD1x,algo::spls->mD1zs);
    //    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodePos-1,1);
    algo::spls->mNodes[mNodePos]->stampV(1,mIndiceStartLambda+1,algo::spls->mD1x,algo::spls->mD1zs);
  }
  if (mNodeNeg >0){
    //    algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodeNeg-1,-1);
    algo::spls->mNodes[mNodeNeg]->stampV(-1,mIndiceStartLambda,algo::spls->mD1x,algo::spls->mD1zs);
    //    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeNeg-1,-1);
    algo::spls->mNodes[mNodeNeg]->stampV(-1,mIndiceStartLambda+1,algo::spls->mD1x,algo::spls->mD1zs);
  }
  //I*lambda
  algo::spls->mD1l->setValue(mIndiceStartLambda,mIndiceStartLambda,1);
  algo::spls->mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda+1,1);
  //+-epsilon
  algo::spls->mD1s->setValue(mIndiceStartLambda,mEpsilon);
  algo::spls->mD1s->setValue(mIndiceStartLambda+1,-mEpsilon);

  //Zns = Vplus +(d11 d12)lambda
  algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,mD11);
  algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+1,mD12);
  algo::spls->mC1s->setValue(mIndiceStartZns,mV2);
}
void componentCOMP_SE1::print(){
  component_LINEAR_NS::print();
  printf("NodeS %d epsilon %f V1 %f V2 %f\n",mNodeS,mEpsilon,mV1,mV2);
  
}
// void componentCOMP_SE1::stampMNA_V(){


//   int i=mI->mIndex;
//   //stamp equations.
//   algo::spls->KCL(mNodeS)->mCoefs[i]=1;
//   //because ie = ie'=0

//   //Zns = Vj-Vk mB1..
//   //VD laws
//   i= algo::spls->getIndexUnknown(ACE_TYPE_V,mNodeS);
//   mEquation->mCoefs[i]+=+1;
//   i= algo::spls->getIndexUnknown(ACE_TYPE_V,0);
//   mEquation->mCoefs[i]-=1;
//   mEquation->mCoefs[mVns->mIndex]-=1;

//   //Y=Vp-Vn + I*lambda +- epsilon
//   //Vp-Vn
//   if (mNodePos >0){
//     algo::spls->mD1x->setValue(mIndiceStartLambda,mNodePos-1,1);
//     algo::spls->mD1x->setValue(mIndiceStartLambda+1,mNodePos-1,1);
//   }
//   if (mNodeNeg >0){
//     algo::spls->mD1x->setValue(mIndiceStartLambda,mNodeNeg-1,-1);
//     algo::spls->mD1x->setValue(mIndiceStartLambda+1,mNodeNeg-1,-1);
//   }
//   //I*lambda
//   algo::spls->mD1l->setValue(mIndiceStartLambda,mIndiceStartLambda,1);
//   algo::spls->mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda+1,1);
//   //+-epsilon
//   algo::spls->mD1s->setValue(mIndiceStartLambda,mEpsilon);
//   algo::spls->mD1s->setValue(mIndiceStartLambda+1,-mEpsilon);

//   //Zns = Vplus +(d11 d12)lambda
//   algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,mD11);
//   algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+1,mD12);
//   algo::spls->mC1s->setValue(mIndiceStartZns,mV2);

  
// }
