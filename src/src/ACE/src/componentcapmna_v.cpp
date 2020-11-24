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

/*! \file componentcapmna_v.cpp

*/
/************************************************************************
  			componentcap.cpp - 


**************************************************************************/

/*



   n       ||      p
---|--->---||------|------
       i   ||

i = Cd(Vn-Vp)/dt is used for the KCL


*/



#include "componentcapmna_v.h"
#include "algo.h"

componentCAPMNA_V::componentCAPMNA_V(dataCAP *d)
:componentCAP(d){
}
componentCAPMNA_V::~componentCAPMNA_V(){
}


void componentCAPMNA_V::stampBeforeInvertion(){
  stamp();
}
void componentCAPMNA_V::stamp(){
  
  //KCL equations
  //  algo::spls->KCL(mData.nodePos)->mCoefs[mU->mDynIndex]-=mData.value;
  algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]-=mData.value;
  algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodePos)]+=mData.value;
 
  //algo::spls->KCL(mData.nodeNeg)->mCoefs[mU->mDynIndex]+=mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->getDynIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=mData.value;

  //U=Vi-Vj
  /*mTenEq->mCoefs[mU->mIndex]-=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg)]+=1;
  mTenEq->mCoefs[algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos)]-=1;*/
}
void componentCAPMNA_V::addUnknowns(){
  //  addTensionUnknown();
}
void componentCAPMNA_V::addEquations(){
  //  addTensionEquation();
}
