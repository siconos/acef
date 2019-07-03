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

/*! \file componentdio.cpp

*/
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
:component_LINEAR_NS(){
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
  mThreshold=ACE_DIODE_THRESHOLD;
  cout<<"*********************************************************************"<<endl;
  cout<<"*************************WARNING, DIODE THRESHOLD IS : "<<mThreshold<<endl;
  

  
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
  ACE_CHECK_IERROR(mIndiceStartLambda == mIndiceStartZns,"dio model component, mIndiceStartLambda != mIndiceStartZns");
  //Bug since the begining, previous code : A lot of test has been validated is the following wrong code :
  //   if (mNodeNeg)
  //     algo::spls->mD1zs->setValue(mIndiceStartZns,mNodeNeg-1,1);
  //   if (mNodePos)
  //     algo::spls->mD1zs->setValue(mIndiceStartZns,mNodePos-1,-1);

  //may be mIndiceStartLambda == mIndiceStartZns. Else big problem.
  if (ACE_FORMULATION !=  ACE_FORMULATION_MNA_V && ACE_FORMULATION !=  ACE_FORMULATION_STAMP_ONLY){
    if (mNodeNeg)
      algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodeNeg-1,1);
    if (mNodePos)
      algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodePos-1,-1);
  }else{
    if (mNodeNeg)
      algo::spls->mD1x->setValue(mIndiceStartLambda,mNodeNeg-1,1);
    if (mNodePos)
      algo::spls->mD1x->setValue(mIndiceStartLambda,mNodePos-1,-1);
  }
  algo::spls->mD1s->setValue(mIndiceStartLambda,mThreshold);
  
}

componentDIO::~componentDIO(){
  
}
void componentDIO::print (){
  component_LINEAR_NS::print();
   printf("diode threshold : %lf \n",mThreshold);
  ;
}
