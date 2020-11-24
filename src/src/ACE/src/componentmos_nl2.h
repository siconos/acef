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

/*! \file componentmos_nl2.h

*/
/**
MOS component:\n
Spice syntax: \n
.model mosP_Sah PMOS LEVEL=1 KP=3.24e-5 VT0=0.6\n
mosP 9 7 0 0 mosP_Sah \n


A model transistor is a necessary information. It specifies witch kind of transistor is used.\n
LEVEL=1 is a necessary option for the parser, but is not used in ACEF.\n
VT0 is the zero-bias voltage.\n
KP is the transconductance parameter.\n
The others parameters of the model are ignored.\n
ACEF use a piecewise linear model for the transistor. The number of hyperplan is a parameter available in the main program.\n

ace formulation:
A Non-Linear smooth formulation.


*/


/************************************************************************
  			componentmos_nl2.h 
**************************************************************************/

#ifndef COMPONENTMOS_NL2_H
#define COMPONENTMOS_NL2_H
#include "ace.h"

#include "component_none_linear_ns.h"

// Class componentMOS
// 
// 
class componentMOS_NL2 : public component_NONE_LINEAR_NS {
 dataMOS1 mData;
  int mNodeD;
  int mNodeG;
  int mNodeS;
  int mMode;
  ACE_DOUBLE mB;
  ACE_DOUBLE mK;
  ACE_DOUBLE mVt;
  equation_NL *mILaw;
public:
  componentMOS_NL2(dataMOS1 *d);
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void  stamp ();
  void stampMNA_V();
  virtual ~componentMOS_NL2();
  virtual void print();

  virtual void computeNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosVector& SICONOS_H);
  virtual void computeJacNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosMatrix& SICONOS_D);


protected:
private:
};
#endif //componentMOS_NL2

