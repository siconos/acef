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

/*! \file componentvsrc_SE1.h

*/
/************************************************************************
  			componentvsrc_SE1.h 

**************************************************************************/

#ifndef COMPONENTVSRC_SE1_H
#define COMPONENTVSRC_SE1_H
#include "componentlinear.h"
#include "ace.h"
// Class componentVSRC
// 
// 
class componentVSRC_SE1 : public componentLINEAR {

public:
  dataVSRC mData;
  ACE_DOUBLE mCurrentValue;
  componentVSRC_SE1(dataVSRC *d);
  virtual ~componentVSRC_SE1();
  virtual void  stamp ();
  virtual void stampTimer();
  void stampTime();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  double * mpCurValue;
protected:
private:
};
#endif //COMPONENTVSRC_SE1_H

