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

/*! \file componentarb.h

*/
/************************************************************************
  			componentarb.h 

**************************************************************************/

#ifndef COMPONENTARB_H
#define COMPONENTARB_H
#include "componentlinear.h"
#include "ace.h"
// Class componentARB
// 
// 
class componentARB : public componentLINEAR {

public:
  dataARB mData;
  ACE_DOUBLE mCurrentValue;
  componentARB(dataARB *d);
  virtual ~componentARB();
  virtual void  stamp ();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
protected:
private:
};
#endif //COMPONENTARB_H

