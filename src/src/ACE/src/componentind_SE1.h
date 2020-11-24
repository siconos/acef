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

/*! \file componentind_SE1.h

*/
/************************************************************************
  			componentind_SE1.h 
**************************************************************************/

#ifndef COMPONENTIND_SE1_H
#define COMPONENTIND_SE1_H
#include "ace.h"
#include "componentdyn.h"

// Class componentIND_SE1
// 
// 
class componentIND_SE1 : public componentDYN {

public:
  dataIND mData;
  componentIND_SE1(dataIND * d);

  virtual void addUnknowns();
  virtual void addEquations();
  virtual void stamp();

virtual ~componentIND_SE1();

protected:
  
private:
  
};
#endif //COMPONENTIND_SE1_H

