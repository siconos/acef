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

/*! \file componentisrc.cpp

*/
/************************************************************************
  			componentisrc.cpp 
**************************************************************************/

#include "componentisrc.h"
#include "algo.h"

componentISRC::~componentISRC(){;}

componentISRC::componentISRC(dataISRC *d)
:componentLINEAR(){
  if(!d)
    ACE_ERROR("ISRC no data");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  mCurrentValue=0;

  mType = ACE_TYPE_ISRC;
  if (ACE_IS_NULL(mData.value))
    ACE_WARNING("ISRC null");
}
void componentISRC::stamp(){
  //nothing, time dependant.
  //KCL
  //  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->mRS]-=mData.value;
  //algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->mRS]+=mData.value;
}


void componentISRC::stampTimer(){
  ACE_DOUBLE newValue;
  ACE_times[ACE_TIMER_TEST_7].start();
  // ParserGetSourceValue("Isource",mData.id,&newValue);
  getISRCValue(mData.id,&newValue);

  //KCL
  algo::spls->KCL(mData.nodeNeg)->mCoefs[algo::spls->mRS]-=newValue - mCurrentValue;
  algo::spls->KCL(mData.nodePos)->mCoefs[algo::spls->mRS]+=newValue - mCurrentValue;
  mCurrentValue = newValue;
  ACE_times[ACE_TIMER_TEST_7].stop();
}



void componentISRC::print(){
  componentLINEAR::print();
  printf("\t value: %f\n",mData.value);
}
