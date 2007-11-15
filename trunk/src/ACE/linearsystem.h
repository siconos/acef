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
  equationVD* addVdEquation();
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
  void set2Sources();

  //DISCRETISATION simulation
  void readInitialValue();
  void initSimu();
  void preparStep();
  bool step();
  void stopSimu();
  void computeZnstiFromX_Zs();
  void simulate();
    
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
  aceMatrix *ms;

  //x'=A1x*x + A1zs*Zs + A1zns*Zns + A1s
  aceMatrix *mA1x;
  aceMatrix *mA1zs;
  aceMatrix *mA1zns;
  aceMatrix *mA1s;
  //0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
  aceMatrix *mB1x;
  aceMatrix *mB1zs;
  aceMatrix *mB1zns;
  aceMatrix *mB1s;

  //Zns = C1x*x + C1s*Zs + C1l*lamdba + C1s
  aceMatrix *mC1x;
  aceMatrix *mC1zs;
  aceMatrix *mC1l;
  aceMatrix *mC1s;

  //Y = D1x*x + D1s*Zs + D1ns*Zns + D1l*lambda +D1s
  aceMatrix *mD1x;
  aceMatrix *mD1zs;
  aceMatrix *mD1zns;
  aceMatrix *mD1l;
  aceMatrix *mD1s;

  aceMatrix *mR;
  aceMatrix *mA2x;
  aceMatrix *mA2zs;
  aceMatrix *mA2s;
  
  aceMatrix *mB2x;
  aceMatrix *mB2zs;
  aceMatrix *mB2l;
  aceMatrix *mB2s;

  aceMatrix *mD2x;
  aceMatrix *mD2zs;
  aceMatrix *mD2l;
  aceMatrix *mD2s;

  

  
  //DISCRETISATION
  aceMatrix *mxti;
  aceMatrix *mzsti;
  aceMatrix *mznsti;
  aceMatrix *mxfree;

  aceMatrix *mW;
  aceMatrix *mD3l;
  aceMatrix *mD3zs;
  aceMatrix *mB3l;
  aceMatrix *mB3zs;
  aceMatrix *mPfree;
  aceMatrix *mQfree;

  ACE_DOUBLE mTheta;
  ACE_DOUBLE mThetap;
  ACE_DOUBLE mH;

  mlcp* mMLCP;
  int mStepCmp;
  int mStepNumber;
  
  ofstream* mSimuStream;
  char mSimuFile[ACE_CHAR_LENGTH];
  
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

protected:
private:
  void buildABCDs();
  void extractDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd);
  void extractNonDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd);
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

  
};
#endif //LINEARSYSTEM_H

