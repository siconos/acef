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

/*! \file algo.h

*/
/************************************************************************
  			algo.h 

**************************************************************************/
/**The algo class.
 */

/**
*/
#ifndef ALGO_H
#define ALGO_H
#include "ace.h"

// Class algo
// 
//
class linearSystem;
class algo {
public:
  static linearSystem *spls;
  static algo *sAlgo;
  /*
   * file : a netlist.
   */
  algo(char * file);
  virtual ~algo();
  void perform();
  void printComponents();
  void stamp();
  void stampAfterInvertion();
  void preparStep(double time);
  void simulate();
  void parseComponents();
  void parseComponents_SE();
  
  void computeNonLinearEquations(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosVector& SICONOS_H); 
  void computeNonLinearJacL_H(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosMatrix& SICONOS_D);




  components mInds;
  components mCaps;
  components mRess;
  components mIsrcs;
  components mVsrcs;
  components mVcvs;
  components mVccs;
  components mArbs;
  components mDios;
  components mMos;
  component_NLs mMos_NL;
  components mBjt;
  components mComps;
  components mRelays;
  char mFile[ACE_CHAR_LENGTH];
  ofstream* mSimuStream;
  char mSimuFile[ACE_CHAR_LENGTH];

protected:
  void performSemiExplicit();
  void performWithOutInvert();
  void performMNA();
  void performMNA_V();
  void performSE();
private:
    

};
#endif //ALGO_H

