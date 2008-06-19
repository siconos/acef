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

  virtual bool step();
  virtual void buildMLCP();
  virtual void fillMLCP();

  void preparMLCP();


protected:
  
};
#endif //LINEARSYSTEMDAE_H

