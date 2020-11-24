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

/*! \file componentdio.h

*/
/**

Diode component.\n
Spice syntax : D1 np nn DIODEF\n

acef formulation:


\f[
Zns=(Id)
\f]
\f[
Zns=\lambda
\f]

\f[
Y=Vn-Vp
\f]

\f[
0 \leq Y \, \perp \, \lambda \geq 0
\f]

*/
/************************************************************************
  			componentdio.h 
**************************************************************************/

#ifndef COMPONENTDIO_H
#define COMPONENTDIO_H
#include "ace.h"

#include "component_linear_ns.h"

// Class componentDIO
// 
// 
class componentDIO : public component_LINEAR_NS {
 dataDIO mData;
  ACE_DOUBLE mThreshold;
public:
  componentDIO(dataDIO *d);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  virtual ~componentDIO();
  virtual void print ();
protected:
private:
};
#endif //COMPONENTDIO_H

