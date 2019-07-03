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

/*! \file componentvcvs_SE1.h

*/
/************************************************************************
  			componentvcvs_SE1.h 

**************************************************************************/

#ifndef COMPONENTVCVS_SE1_H
#define COMPONENTVCVS_SE1_H
#include "componentlinear.h"
#include "ace.h"
// Class componentVCVS_SE1
// 
// 
class componentVCVS_SE1 : public componentLINEAR {

public:
  dataVCVS mData;
  ACE_DOUBLE mCurrentValue;
  componentVCVS_SE1(dataVCVS *d);
  virtual ~componentVCVS_SE1();
  virtual void  stamp ();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void print();
protected:
private:
};
#endif //COMPONENTVCVS_SE1_H

