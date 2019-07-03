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

/*! \file componentres.cpp

*/
/************************************************************************
  			componentres.cpp 
**************************************************************************/

#include "componentres.h"
#include "algo.h"

componentRES::componentRES(dataRES *d)
:componentLINEAR(){
  if (!d)
    ACE_ERROR("componentRES::componentRES dataRES null");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  if (ACE_IS_NULL(mData.value))
    ACE_ERROR("Resistor null");
  
  mType = ACE_TYPE_RES;
}
/*stamp:
  
     U=Ri
          U=Vn-Vp
     <-----------
----|------->----|--->
    n       i    p

        | Vp   |  Vn   |
___________________________
KCL(p)  |-1/R  | 1/R   |
___________________________
KCL(n)  |1/R  | -1/R   |
*/
void componentRES::stamp(){
  //KCLs
  int vp,vn;
  //equationKCL *pk;
  vp=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos);
  vn=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  //pk=;
  algo::spls->KCL(mData.nodePos)->mCoefs[vp]+=-1/mData.value;
  algo::spls->KCL(mData.nodePos)->mCoefs[vn]+=1/mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[vp]+=1/mData.value;
  algo::spls->KCL(mData.nodeNeg)->mCoefs[vn]+=-1/mData.value;

}
void componentRES::print(){
  componentLINEAR::print();
  printf("\t value: %f\n",mData.value);
}

componentRES::~componentRES(){
  
}
