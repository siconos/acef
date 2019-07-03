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

/*! \file equation.cpp

*/
/************************************************************************
  			equation.cpp 
**************************************************************************/
#include "equation.h"
#include <iostream>
using namespace std;

equation::equation(){
  mLine=-1;
  mCoefs=(ACE_DOUBLE *)0;
  mIsDyn = false;
  mIsLinear = true;
  mAvailable=true;
  mSize = 0;

}

void equation::allocMemory(int nb)
 {
   if(nb<1){
     ACE_INTERNAL_ERROR("equation alloc with nb<1");
     return;
   }
   mSize = nb;
    mCoefs = (ACE_DOUBLE*)calloc(mSize,sizeof(ACE_DOUBLE));
 }

 void equation::print(){
   for(int i =0; i <mSize;i++)
     cout <<"\t"<<mCoefs[i];
   cout <<"\n";
 }

equation::~equation(){
  if (mCoefs)
    free(mCoefs);
}


