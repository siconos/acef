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

/*! \file componentcapdae.cpp

*/
/************************************************************************
  			componentcap.cpp - 


**************************************************************************/

/*



   n       ||      p
---|--->---||------|------
       i   ||

Unp = Vn-Vp

Cd(Unp)/dt = i

i is added in the vector of unknowns. i is used to fill the KCL

*/



#include "componentcapdae.h"
#include "algo.h"

componentCAPDAE::componentCAPDAE(dataCAP *d)
:componentCAP(d){
}
componentCAPDAE::~componentCAPDAE(){
}


void componentCAPDAE::stampBeforeInvertion(){
  stamp();
}
void componentCAPDAE::stamp(){
  ACE_CHECK_IERROR(mU,"componentCAPDAE::stamp mU null");
  ACE_CHECK_IERROR(mI,"componentCAPDAE::stamp mI null");
  ACE_CHECK_IERROR(mTenEq,"componentCAPDAE::stamp mTenEq null");
  ACE_CHECK_IERROR(mDynEquation,"componentCAPDAE::stamp mDynEquation null");
  
  bool dyn=false;
  //KCL equations
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]-=1;
  algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;

  //U=Vi-Vj
  mTenEq->mCoefs[mU->mIndex]-=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;
  
  //Cu'=I
  mDynEquation->mCoefs[mU->mDynIndex]+=mData.value;
  mDynEquation->mCoefs[mI->mIndex]+=1;
  
}
