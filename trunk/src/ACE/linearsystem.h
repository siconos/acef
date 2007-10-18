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

// Class linearSystem
// 
// 
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
  void print();

  unknown* addinx(int type, component* c);
  unknown* addinZs(int type, component* c);
  unknown* addinZns(int type, component* c);
  equationKCL* KCL(int i);
  equationCAP* addCapEquation();
  equationIND* addIndEquation();
  equationVD* addVdEquation();
  equationTEN* addTenEquation();
  int getIndexUnknown (int type,int node);

  int mNbNodes;
  int mNbUnknowns;
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
protected:
private:
};
#endif //LINEARSYSTEM_H

