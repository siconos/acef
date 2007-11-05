/************************************************************************
  			equation.h 
**************************************************************************/

#ifndef EQUATION_H
#define EQUATION_H
#include "ace.h"

class equation {
public:
  equation();
  virtual ~equation();
  void allocMemory(int nb);
  virtual void print();
  ACE_DOUBLE* mCoefs;
  int mLine;
  bool mIsDyn;
  bool mAvailable;
  int mSize;

protected:
private:
};
#endif //EQUATION_H

