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


  int mNbNodes;
  int mNbUnknowns;
  int mNbEquations;
  int mNbDynEquations;
  int mNbNonDynEquations;

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


  
  void printEquations(ostream& os = cout);
  void printABCDs(ostream& os = cout);
  void printA1(ostream& os = cout);
  void printB1(ostream& os = cout);
  void printSystemInTabFile(char * file);

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

};
#endif //LINEARSYSTEM_H

