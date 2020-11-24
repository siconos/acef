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

/*! \file componentcomp_SE1.h

*/
/**

Comparator component.\n
ACEF syntax : .comp N+ N- Noutput Vmin=0 Vmax=3 Vepsilon=0.1\n

characteristic:\n


\image html comparator.jpg


acef formulation:\n
\f[
dim(\lambda) = 2
\f]
\f[
Zns=Vs
\f]
\f[
Y = \left(\begin{array}{c}
Y_1\\
Y_2\\
\end{array}\right)=

\left(\begin{array}{c}
Vp-Vn\\
Vp-Vn\\
\end{array}\right)+I\lambda+
\left(\begin{array}{c}
epsilon\\
-epsilon\\
\end{array}\right)
\f]
\f[
Zns = Vplus + (d11,d12)\lambda
\f]
*/

/************************************************************************
  			componentcomp_SE1.h 

**************************************************************************/

#ifndef COMPONENTCOMP_SE1_H
#define COMPONENTCOMP_SE1_H
#include "ace.h"
#include "component_linear_ns.h"

// Class componentCOMP
// 
// 
class componentCOMP_SE1 : public component_LINEAR_NS {
public:
  dataCOMP mData;
  ACE_DOUBLE mV2;
  ACE_DOUBLE mV1;
  ACE_DOUBLE mEpsilon;
  ACE_DOUBLE mD11;
  ACE_DOUBLE mD12;
  int mNodeS;
  unknown *mVns;
  
  componentCOMP_SE1(dataCOMP *d);
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void  stamp ();
  //  void stampMNA_V();

  virtual void  print ();

  virtual ~componentCOMP_SE1();
protected:
  
private:
 
};
#endif //COMPONENTCOMP_SE1_H

