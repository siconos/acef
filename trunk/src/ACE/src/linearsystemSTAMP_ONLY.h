/************************************************************************
  			linearsystemSTAMP_ONLY.h 
**************************************************************************/

#ifndef LINEARSYSTEMSTAMP_ONLY_H
#define LINEARSYSTEMSTAMP_ONLY_H
#include "linearsystemMNA_V.h"
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
class linearSystemSTAMP_ONLY : public linearSystemMNA_V {
public:
  linearSystemSTAMP_ONLY();
  virtual ~linearSystemSTAMP_ONLY();
  virtual void buildMLCP();
  virtual void fillMLCP();
  virtual void preparMLCP();
  virtual bool step();
  virtual void extractDynamicSystemSource();
  virtual void ExtractAndCompute2Sources();

private:
  aceMatrix *mA_A1;
  aceVector *mA1sti;

};
#endif //LINEARSYSTEMSTAMP_ONLY_H

