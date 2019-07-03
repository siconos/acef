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

/*! \file componentcapmna.cpp

*/
/************************************************************************
  			componentcap.cpp - 


**************************************************************************/

/*



   n       ||      p
---|--->---||------|------
       i   ||

Unp = Vn-Vp

i = Cd(Unp)/dt is used for the KCL


*/



#include "componentcapmna.h"
#include "algo.h"

componentCAPMNA::componentCAPMNA(dataCAP *d)
:componentCAP(d){
}
componentCAPMNA::~componentCAPMNA(){
}


void componentCAPMNA::stampBeforeInvertion(){
  stamp();
}
void componentCAPMNA::stamp(){
  ACE_CHECK_IERROR(mU,"componentCAPMNA::stamp mU null");
  ACE_CHECK_IERROR(mTenEq,"componentCAPMNA::stamp mTenEq null");
  
  bool dyn=false;
  //KCL equations
  algo::spls->KCL(mData.nodePos)->mCoefs[mU->mDynIndex]-=mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mU->mDynIndex]+=mData.value;
  
  //U=Vi-Vj
  mTenEq->mCoefs[mU->mIndex]-=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;
}
void componentCAPMNA::addUnknowns(){
  addTensionUnknown();
}
void componentCAPMNA::addEquations(){
  addTensionEquation();
}
