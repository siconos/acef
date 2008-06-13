/************************************************************************
  			linearsystemDAE.h 
**************************************************************************/

#ifndef LINEARSYSTEMDAE_H
#define LINEARSYSTEMDAE_H
#include "linearsystem.h"
using namespace std;
// Class linearSystem
// 
// Ax'=Bx+CZs+DZns+s
class linearSystemDAE : public linearSystem {
public:
  linearSystemDAE();
  virtual ~linearSystemDAE();

  virtual void preparStep();
  virtual bool step();
  virtual void buildMLCP();
  
};
#endif //LINEARSYSTEMDAE_H

