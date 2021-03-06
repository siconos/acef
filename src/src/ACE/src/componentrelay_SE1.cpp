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

/*! \file componentrelay_SE1.cpp

*/
/************************************************************************
  			componentrelay_SE1.cpp 
**************************************************************************/

#include "componentrelay_SE1.h"
#include "algo.h"

componentRELAY_SE1::componentRELAY_SE1(dataCOMP *d)
:component_LINEAR_NS(){
  ACE_CHECK_IERROR(d,"componentRELAY_SE1::componentRELAY_SE1 : Diode data null");
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
componentRELAY_SE1::~componentRELAY_SE1(){;}

void componentRELAY_SE1::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
  mVns=algo::spls->addinZns(ACE_TYPE_U,this);
  mIndiceStartZns= mVns->mIndexInVector;
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}

void componentRELAY_SE1::addEquations(){
  mEquation=algo::spls->addVdEquation(mName);
}
void componentRELAY_SE1::stamp(){
  int coef=1;
  if (mV2>mV1)
    coef=-1;
  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mNodeS)->mCoefs[i]=1;
  //because ie = ie'=0

  //VD laws
  //Zns = Vs-V0 ..
  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mNodeS);
  //  mEquation->mCoefs[i]+=+1;
  algo::spls->mNodes[mNodeS]->stampV(1,mEquation->mCoefs + algo::spls->mDimx,mEquation->mCoefs + 2*algo::spls->mDimx);
  
  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,0);
  //  mEquation->mCoefs[i]-=1;
  algo::spls->mNodes[0]->stampV(-1,mEquation->mCoefs + algo::spls->mDimx,mEquation->mCoefs + 2*algo::spls->mDimx);
  
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
    //    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodePos-1,-1);
    algo::spls->mNodes[mNodePos]->stampV(-1,mIndiceStartLambda+1,algo::spls->mD1x,algo::spls->mD1zs);
  }
  if (mNodeNeg >0){
    //    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeNeg-1,1);
    algo::spls->mNodes[mNodeNeg]->stampV(1,mIndiceStartLambda+1,algo::spls->mD1x,algo::spls->mD1zs);
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
void componentRELAY_SE1::print(){
  component_LINEAR_NS::print();
  printf("NodeS %d offset %f V1 %f V2 %f\n",mNodeS,mOffset,mV1,mV2);
  
}
