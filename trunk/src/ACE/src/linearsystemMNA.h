/************************************************************************
  			linearsystemMNA.h 
**************************************************************************/

#ifndef LINEARSYSTEMMNA_H
#define LINEARSYSTEMMNA_H
#include "linearsystem.h"
using namespace std;
// Class linearSystem
//
// MX'= AX + R*lambda +s
// Y = DX + E*lambda +s
// 0<Y ortho lambda >0
//
// Ax'=Bx+CZs+DZns+s
class linearSystemMNA : public linearSystem {
public:
  linearSystemMNA();
  virtual ~linearSystemMNA();

  virtual bool step();
  virtual void buildMLCP();
  virtual void fillMLCP();
  virtual void extractDynamicSystemSource();

  
  void preparMLCP();

  virtual void printA1(ostream& os );
  virtual void printSystem2(ostream& os = cout);

protected:
  virtual void allocMemory();
};
#endif //LINEARSYSTEMMNA_H

