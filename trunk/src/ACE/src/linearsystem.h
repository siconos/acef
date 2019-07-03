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

/*! \file linearsystem.h

*/
/************************************************************************
  			linearsystem.h 
**************************************************************************/

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H
#include "unknown.h"
#include "ace.h"
#include "equationkcl.h"
#include "equationcap.h"
#include "equationind.h"
#include "equationvd.h"
#include "equationten.h"
#include "equation_nl.h"
#include "aceMatrix.h"
#include "aceVector.h"
#include "mlcp.h"
#include "node.h"
using namespace std;
// Class linearSystem
// 
// Ax'=Bx+CZs+DZns+s
class linearSystem {
public:
  linearSystem();
  virtual ~linearSystem();
  void preparForStamp();
  virtual void allocMemory();
  void initKCL();
  virtual void addVUnknowns();
  bool isUnknown (int type, component* c,unknown **uout);
  void addKCLinDyn(int j);

  nodes mNodes;

  unknown* addinx(int type, component* c);
  unknown* addinZs(int type, component* c);
  unknown* addinZns(int type, component* c);
  equationKCL* KCL(int i);
  equationCAP* addCapEquation();
  equationIND* addIndEquation();
  equationVD* addVdEquation(char* name =0);
  equationTEN* addTenEquation();
  equation_NL* addNonLinearEquation();
  virtual int getIndexUnknown (int type,int node);
  virtual int getDynIndexUnknown (int type,int node){;};

  //write x'=A1x * mx + A1zs * mZs + A1zns * mZns;
  void computedxdt();
  void getlinefromdxdt(int line, ACE_DOUBLE * coefs);

  //0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
  void buildLinearSystem();

  void buildABCDs();

  //
  //R=A1zns*C1l
  //x'=A2x*x + A2zs*Zs + R*lambda+A2s
  //0=B2x*x + B2zs*Zs + B2l*lambda + B2s
  //Y=D2x*x + D2zs*Zs + D2l*lambda + D2s
  //
  virtual void alloc2Matrix();
  virtual void set2matrix();

  //DISCRETISATION simulation
  void readInitialValue();
  virtual double getCurrentTime();
  virtual void initSimu(ofstream* pstream);
  virtual void preparStep();
  /*compute the right side of the mlcp*/
  virtual void preparMLCP();
  virtual bool step();
  void stopSimu();
  void printLog();
  virtual void computeZnstiFromX_Zs();
  /*build the mlcp:
  *
  *IF ADAPTIVE TIME STEPPING
  * If PRE_COMPUTE_MLCP
  *  call only one time, and precompute system for each possible step (see ace.h)
  * ELSE
  *  call at each step when the step change.
  * END
  *ELSE
  * call only on time
  *END
  */
  virtual void buildMLCP();
  /*
   *Fill the MLCP matrix with current step
   *
   */
  virtual void fillMLCP();
  /*compute the largest time stepping that respect the ACE_MAX_LOCAL_ERROR*/
  virtual void computeAndAcceptStep();
  virtual void computeAndAcceptStep_2();

  /**
     These fucntions are called at each step to update de terme source (Voltage source and Current source)
   */
  virtual void ExtractAndCompute2Sources();
  void extractSources();

  virtual void extractDynamicSystemSource();
  virtual void extractInteractionSource();
  virtual void updateDynamicSystemSource();
  virtual void updateInteractionSource();
  
  long  mLogFrequency;
  long  mLogPrint;
  long mPourMille;

    
  char mFile[ACE_CHAR_LENGTH];

  int mNbNodes;
  int mNbUnknowns;
  int mNbEquations;
  int mNbDynEquations;
  int mNbNonDynEquations;

  int mDimLambda;
  int mDimx;
  int mDimzs;
  int mDimzns;
  /*mV0x = 1 if V0 is in x else mV0x=0*/
  int mV0x;
  /*mV0zs = 1 if V0 is in zs else mV0zs=0*/
  int mV0zs;

  int mRS;
  bool mReAlloc;
  int mDynIndex;
  unknowns mx;
  unknowns mZs;
  unknowns mZns;
  equations mKCL;
  equations mVD;
  equations mIND;
  equations mCAP;
  equations mTEN;
  equations mNL;

  //siconos obj
  // Ax'=Bx+CZs+DZns+s
  aceMatrix *mA;
  aceMatrix *mB;
  aceMatrix *mC;
  aceMatrix *mD;
  aceVector *ms;
  aceMatrix *mMatBuf1;
  aceMatrix *mMatBuf2;

  //x'=A1x*x + A1zs*Zs + A1zns*Zns + A1s
  aceMatrix *mA1x;
  aceMatrix *mA1zs;
  aceMatrix *mA1zns;
  aceVector *mA1s;
  //0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
  aceMatrix *mB1x;
  aceMatrix *mB1zs;
  aceMatrix *mB1zns;
  aceVector *mB1s;

  //Zns = C1x*x + C1s*Zs + C1l*lamdba + C1s
  aceMatrix *mC1x;
  aceMatrix *mC1zs;
  aceMatrix *mC1l;
  aceVector *mC1s;

  //Y = D1x*x + D1s*Zs + D1ns*Zns + D1l*lambda +D1s
  aceMatrix *mD1x;
  aceMatrix *mD1zs;
  aceMatrix *mD1zns;
  aceMatrix *mD1l;
  aceVector *mD1s;

  aceMatrix *mR;
  aceMatrix *mhR[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mA2x;
  aceMatrix *mA2zs;
  aceMatrix *mHThetaA2zs[ACE_NB_ADAPT_STEP+1];
  
  aceVector *mA2s;
  aceVector *mA2sti;

  aceMatrix *mB2x;
  aceMatrix *mB2zs;
  aceMatrix *mB2l;
  aceVector *mB2s;

  aceMatrix *mD2x;
  aceMatrix *mD2zs;
  aceMatrix *mD2l;
  aceVector *mD2s;

  aceMatrix *mD2xW[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mB2xW[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mHThetaWA2zs[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mHWR[ACE_NB_ADAPT_STEP+1];


  //Simulation
  ofstream* mSimuStream;

  
  //DISCRETISATION
  aceVector *mxti;
  aceVector *mzsti;
  aceVector *mznsti;
  aceVector *mlti;
  aceVector *mxfree;

 

  
  aceMatrix *mW[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mD3l[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mD3zs[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mB3l[ACE_NB_ADAPT_STEP+1];
  aceMatrix *mB3zs[ACE_NB_ADAPT_STEP+1];
  aceVector *mPfree;
  aceVector *mQfree;
  aceVector *mPAux;
  aceMatrix *mPxAux;
 

  ACE_DOUBLE mTheta;
  ACE_DOUBLE mThetazs;
  ACE_DOUBLE mThetap;
  ACE_DOUBLE mH;
  ACE_DOUBLE mHori;
  ACE_DOUBLE mMaxHori;
  //for adaptive time stepping
  bool mUseAdaptiveTimeStepping;
  ACE_DOUBLE mAlpha;
  ACE_DOUBLE mAlphaMax;
  ACE_DOUBLE mAlphaMin;
  ACE_DOUBLE mNormX0;
  ACE_DOUBLE mNormZ0;
  int mAdaptCmp;
  int mAllStepCmp;

  ACE_DOUBLE mTstart;
  ACE_DOUBLE mTstop;
  ACE_DOUBLE mTcurrent;
  aceVector *mxtiprev;
  aceVector *mzstiprev;
  aceVector *mxticurrent;
  aceVector *mzsticurrent;
  aceVector *mzst_inter;
  aceVector *mxt_inter;
  aceVector *mzst1;
  aceVector *mxt1;
  aceVector *mxbuf;
  aceVector *mzsbuf;
  

  mlcp* mMLCP;
  long mStepCmp;
  long mStepNumber;
  
  
  void printEquations(ostream& os = cout);
  void printABCDs(ostream& os = cout);
  virtual void printA1(ostream& os = cout);
  void printB1(ostream& os = cout);
  void printC1(ostream& os = cout);
  void printD1(ostream& os = cout);
  void printSystemInTabFile(char * file);
  virtual void printSystem2(ostream& os = cout);
  virtual void printStep(ostream& os,aceVector *pVx,aceVector *pVzs);
  ACE_DOUBLE computeAnalyticSolution(ACE_DOUBLE t);
  ACE_DOUBLE mSommeError;
  ACE_DOUBLE mCoef;

  
  void printDiscretisation(ostream& os = cout);

  virtual void allocForInitialValue();


protected:
  virtual void setStep(ACE_DOUBLE newH);
  bool mAdaptiveStepEvaluation;

  ACE_DOUBLE mSerrorX;
  ACE_DOUBLE mSerrorZs;
  ACE_DOUBLE mMaxError;
  int mNbToSmall;
  int mNbToBig;
  int mNbBacktrack;
  void extractDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd);
  void extractDynBockInVect(aceVector * V);


  void extractNonDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd);
  void extractNonDynBockInVect(aceVector * V);

  void allocA1Matrix();
  void freeA1Matrix();
  void allocB1Matrix();
  void freeB1Matrix();
  void allocC1Matrix();
  void freeC1Matrix();
  void allocD1Matrix();
  void freeD1Matrix();
  virtual void allocDiscretisation();
  void freeDiscretisation();
  void freeForInitialValue();
  
  
};
#endif //LINEARSYSTEM_H

