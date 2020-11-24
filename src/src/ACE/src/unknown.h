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

/*! \file unknown.h

*/
/************************************************************************
  			unknown.h 
**************************************************************************/

#ifndef UNKNOWN_H
#define UNKNOWN_H
#include "ace.h"

#include "component.h"
// Class unknown
// 
// 
class unknown {
public:
   unknown(int type, component *c);
   unknown(int type, int node);
  void print();
  void printdev();
  int mNode;
   int mType;
  //mIndex in ----x'-------x--------Zs---------Zns--------
   int mIndex;
  //mDynIndex only for unknown in x: ----x'-------x--------Zs---------Zns--------
   int mDynIndex;
  //mIndexInVector index in x or Zs or Zns.
  int mIndexInVector;
  component * mComponent;
  char mName[64];
protected:
private:
  void buildName();
};
#endif //UNKNOWN_H

