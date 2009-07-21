/************************************************************************
  			linearsystemMNA.h 
**************************************************************************/

#ifndef LINEARSYSTEMMNA_V_H
#define LINEARSYSTEMMNA_V_H
#include "linearsystem.h"
using namespace std;
// Class linearSystem
//
// GOAL: DO NOT ADD CAP TENSION IN X
//
//
// MX'= AX + R*lambda +s
// 0 = DX + E*lambda +s
// 0<Y ortho lambda >0
//
// Ax'=Bx+CZs+DZns+s
class linearSystemMNA_V : public linearSystem {
public:
  linearSystemMNA_V();
  virtual ~linearSystemMNA_V();

  virtual bool step();
  virtual void buildMLCP();
  virtual void fillMLCP();
  virtual void extractDynamicSystemSource();

  
  virtual void preparMLCP();

  virtual void printA1(ostream& os );
  virtual void printSystem2(ostream& os = cout);
  virtual void addVUnknowns();
  virtual int getIndexUnknown (int type,int node);
  virtual int getDynIndexUnknown (int type,int node);
  virtual void printStep(ostream& os,aceVector *pVx,aceVector *pVzs);


protected:
  virtual void allocMemory();
};
#endif //LINEARSYSTEMMNA_V_H

