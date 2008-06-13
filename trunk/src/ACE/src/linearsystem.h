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
#include "aceMatrix.h"
#include "aceVector.h"
#include "mlcp.h"
using namespace std;
// Class linearSystem
// 
// Ax'=Bx+CZs+DZns+s
class linearSystem {
public:
  linearSystem();
  virtual ~linearSystem();
  void preparForStamp();
  void allocMemory();
  void initKCL();
  void addVUnknowns();
  bool isUnknown (int type, component* c,unknown **uout);
  void addKCLinDyn(int j);



  unknown* addinx(int type, component* c);
  unknown* addinZs(int type, component* c);
  unknown* addinZns(int type, component* c);
  equationKCL* KCL(int i);
  equationCAP* addCapEquation();
  equationIND* addIndEquation();
  equationVD* addVdEquation(char* name =0);
  equationTEN* addTenEquation();
  int getIndexUnknown (int type,int node);

  //write x'=A1x * mx + A1zs * mZs + A1zns * mZns;
  void computedxdt();
  void getlinefromdxdt(int line, ACE_DOUBLE * coefs);

  //0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
  void buildLinearSystem();


  //
  //R=A1zns*C1l
  //x'=A2x*x + A2zs*Zs + R*lambda+A2s
  //0=B2x*x + B2zs*Zs + B2l*lambda + B2s
  //Y=D2x*x + D2zs*Zs + D2l*lambda + D2s
  //
  void set2matrix();

  //DISCRETISATION simulation
  void readInitialValue();
  virtual double getCurrentTime();
  virtual void initSimu();
  virtual void preparStep();
  virtual bool step();
  void stopSimu();
  virtual void computeZnstiFromX_Zs();
  virtual void buildMLCP();
  virtual void computeBestStep();

  void ExtractAndCompute2Sources();
  void extractSources();

  void extractDynamicSystemSource();
  void extractInteractionSource();
  void computeDynamicSystemSource();
  void computeInteractionSource();
  
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
  aceMatrix *mhR;
  aceMatrix *mA2x;
  aceMatrix *mA2zs;
  aceMatrix *mHThetaA2zs;
  
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

  aceMatrix *mD2xW;
  aceMatrix *mB2xW;
  aceMatrix *mHThetaWA2zs;
  aceMatrix *mHWR;


  
  //DISCRETISATION
  aceVector *mxti;
  aceVector *mzsti;
  aceVector *mznsti;
  aceVector *mxfree;

 

  
  aceMatrix *mW;
  aceMatrix *mD3l;
  aceMatrix *mD3zs;
  aceMatrix *mB3l;
  aceMatrix *mB3zs;
  aceVector *mPfree;
  aceVector *mQfree;
  aceVector *mPAux;
  aceMatrix *mPxAux;

  ACE_DOUBLE mTheta;
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
  ACE_DOUBLE mLocalErrorTol;
  int mAdaptCmp;

  ACE_DOUBLE mTstart;
  ACE_DOUBLE mTstop;
  ACE_DOUBLE mTcurrent;
  aceVector *mxtiprev;
  aceVector *mzstiprev;
  aceVector *mxticurrent;
  aceVector *mzsticurrent;
  

  mlcp* mMLCP;
  long mStepCmp;
  long mStepNumber;
  
  
  void printEquations(ostream& os = cout);
  void printABCDs(ostream& os = cout);
  void printA1(ostream& os = cout);
  void printB1(ostream& os = cout);
  void printC1(ostream& os = cout);
  void printD1(ostream& os = cout);
  void printSystemInTabFile(char * file);
  void printSystem2(ostream& os = cout);
  void printStep(ostream& os = cout);
  void printDiscretisation(ostream& os = cout);

  void allocForInitialValue();


protected:
  virtual void setStep(ACE_DOUBLE newH);
private:
  void buildABCDs();
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
  void allocDiscretisation();
  void freeDiscretisation();
  void freeForInitialValue();
  
  
};
#endif //LINEARSYSTEM_H

