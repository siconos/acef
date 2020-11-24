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

/*! \file component.cpp

*/
/************************************************************************
  			component.cpp 
**************************************************************************/


#include "unknown.h"
#include "equation.h"

component::~component(){
}
component::component(){
  mU=0;
  mI=0;
  mEquation=0;
  mType=ACE_TYPE_NO;
  mNodePos=-1;
  mNodeNeg=-1;
  mName=0;

}

void component::stampTimer(){
  ;}

void component::stamp () {
  
}
/**
 * 
 */
void component::addUnknowns () {
  
}
/**
 * 
 */
void component::addEquations () {
  
}
void component::print () {
  char name[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,name);
  printf("component type : %s \n",name);
  if (mName)
    printf("\t Name : %s\n",mName);
  printf("\tnodePos, nodeNeg : %d %d\n",mNodePos,mNodeNeg);
  
  
}


