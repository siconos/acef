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

/*! \file acefRelation.h

*/
#ifndef ACEFRELATION_H
#define ACEFRELATION_H

#define WITH_KERNEL_RELATION

#ifdef WITH_KERNEL_RELATION

#include "SiconosKernel.hpp"
#include "ace.h"

class acefRelation : public FirstOrderType2R
{
protected:
  SimpleMatrix* mB;
  SimpleMatrix* mC;
  SimpleMatrix* mD;
  ACE_DOUBLE mLastComputedSource;
public:
  acefRelation();
  virtual ~acefRelation(){};


  virtual void initialize(SP::Interaction inter);
  void initJac(SimpleMatrix* C,SimpleMatrix* D,SimpleMatrix* B);

  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeH(double) ;
	
  /** default function to compute g 
   *  \param double : current time
   */
  virtual void computeG(double) ;

  /** default function to compute jacobianH
   *  \param double : current time
   */
  virtual void computeJacXH(double);
  virtual void computeJacLH(double);
	
  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacXG(double);
  virtual void computeJacLG(double);
	


  double source(double t);

};

TYPEDEF_SPTR(acefRelation);
#endif

#endif
