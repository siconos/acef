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

/*! \file componentvsrc_SE1.cpp

*/
/************************************************************************
  			componentvsrc_SE1.cpp 
**************************************************************************/

#include "componentvsrc_SE1.h"
#include "algo.h"

componentVSRC_SE1::~componentVSRC_SE1(){;}
componentVSRC_SE1::componentVSRC_SE1(dataVSRC *d)
:componentLINEAR(){
  if(!d)
    ACE_ERROR("VSRC no data");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mpCurValue = mData.pCurValue;
  mType = ACE_TYPE_VSRC;
  mCurrentValue=0;
  if (ACE_IS_NULL(mData.value))
    ACE_WARNING("VSRC null");
}
void componentVSRC_SE1::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
}
void componentVSRC_SE1::addEquations(){
  mEquation=algo::spls->addVdEquation();
}
void componentVSRC_SE1::stamp(){
  ACE_DOUBLE newValue;

  //KCL
  int i;
  ACE_CHECK_IERROR(mI && mEquation,"componentVSRC_SE1::stamp no mI or no mEquation!!");
  algo::spls->KCL(mData.nodeNeg)->mCoefs[mI->mIndex]-=1;
  algo::spls->KCL(mData.nodePos)->mCoefs[mI->mIndex]+=1;
  
  //VD laws
  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  //mEquation->mCoefs[i]+=-1;
  algo::spls->mNodes[mData.nodeNeg]->stampV(-1,mEquation->mCoefs + algo::spls->mDimx,mEquation->mCoefs + 2*algo::spls->mDimx);
  //  i= algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos);
  //  mEquation->mCoefs[i]+=1;
  algo::spls->mNodes[mData.nodePos]->stampV(1,mEquation->mCoefs + algo::spls->mDimx,mEquation->mCoefs + 2*algo::spls->mDimx);
}

void componentVSRC_SE1::stampTimer(){
  ACE_DOUBLE newValue;
  //  ParserGetSourceValue("Vsource",mData.id,&newValue);
  //  getVSRCValue(mData.id,&newValue);
  newValue = (*mpCurValue);

  mEquation->mCoefs[algo::spls->mRS]-=newValue - mCurrentValue;
  mCurrentValue = newValue;
}


