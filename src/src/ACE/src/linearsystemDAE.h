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

/*! \file linearsystemDAE.h

*/
/************************************************************************
  			linearsystemDAE.h 
**************************************************************************/

#ifndef LINEARSYSTEMDAE_H
#define LINEARSYSTEMDAE_H
#include "linearsystem.h"
using namespace std;
// Class linearSystem
// 
// Ax'=Bx+CZs+DZns+s
class linearSystemDAE : public linearSystem {
public:
  linearSystemDAE();
  virtual ~linearSystemDAE();

  virtual bool step();
  virtual void buildMLCP();
  virtual void fillMLCP();
  virtual void extractDynamicSystemSource();

  
  void preparMLCP();

  virtual void printA1(ostream& os );
  virtual void printSystem2(ostream& os = cout);

protected:
  
};
#endif //LINEARSYSTEMDAE_H

