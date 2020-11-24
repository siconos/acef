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

/*! \file componentcap.h

*/
/************************************************************************
  			componentcap.h
**************************************************************************/

#ifndef COMPONENTCAP_H
#define COMPONENTCAP_H
#include "ace.h"

#include "componentdyn.h"
#include "equationten.h"

// Class componentCAP
// 
//
/**
 
 It implements the capacitor component.\n
 
 The constitutive equation is:\n
 \f[
 U_{np} = V_n-V_p
 \f]
 \f[
 C*\frac{U_{np}}{dt} = i(t)
 \f]
 
 */
class componentCAP : public componentDYN {
public:
  dataCAP mData;
  componentCAP(dataCAP *d);
  ACE_DOUBLE * mICoefs;
  virtual ~componentCAP();
  /*
   *
   *
   */
  virtual void  stampBeforeInvertion ();
  /*
   *
   *
   */
  void  stampAfterInvertion ();
  /*
   *
   *
   */
  virtual void addUnknowns();
  virtual void addEquations();

  /*
   *
   *
   */
  void addCurrentEquation();
  /*
   *
   *
   */
  void addTensionEquation();
  
  /**
   *Add the constitutive equation.
   *
   */
  void addCurrentUnknown();
  /**
   *Add Unp as an unknown.
   *
   */
  void addTensionUnknown();
  equationTEN* mTenEq;

  void printI();

protected:
 
private:
};
#endif //COMPONENTCAP_H



