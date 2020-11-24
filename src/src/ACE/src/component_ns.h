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

/*! \file component_ns.h

*/
/************************************************************************
  			component_NONE_SMOOTH.h 
**************************************************************************/

#ifndef COMPONENT_NS_H
#define COMPONENT_NS_H
#include "ace.h"
#include "component.h"

// Class component_NS
// 
// 
class component_NS : public component {

public:
  component_NS();
  int mDimlambda;
  int mDimZns;
  int mIndiceStartZns;
  int mIndiceStartLambda;
  
protected:

private:
 
};


#endif //COMPONENT_NS_H

