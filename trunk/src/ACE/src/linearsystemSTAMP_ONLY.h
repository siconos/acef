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

/*! \file linearsystemSTAMP_ONLY.h

*/
/************************************************************************
  			linearsystemSTAMP_ONLY.h 
**************************************************************************/

#ifndef LINEARSYSTEMSTAMP_ONLY_H
#define LINEARSYSTEMSTAMP_ONLY_H
#include "linearsystemMNA_V.h"
using namespace std;
// Class linearSystem
//
// GOAL: DO NOT ADD CAP TENSION IN X
//
//
// MX'= AX + R*lambda +s
// 0 = DX + E*lambda +s
// 0<Y ortho lambda >0
//
// Ax'=Bx+CZs+DZns+s
class linearSystemSTAMP_ONLY : public linearSystemMNA_V {
public:
  linearSystemSTAMP_ONLY();
  virtual ~linearSystemSTAMP_ONLY();
  virtual void buildMLCP();
  virtual void fillMLCP();
  virtual void preparMLCP();
  virtual bool step();
  virtual void extractDynamicSystemSource();
  virtual void ExtractAndCompute2Sources();

private:
  aceMatrix *mA_A1;
  aceVector *mA1sti;

};
#endif //LINEARSYSTEMSTAMP_ONLY_H

