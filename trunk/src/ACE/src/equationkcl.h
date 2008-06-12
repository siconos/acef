/************************************************************************
  			equationkcl.h 
**************************************************************************/

#ifndef EQUATIONKCL_H
#define EQUATIONKCL_H
#include "equation.h"
#include "ace.h"

// Class equationkcl
// 
// 
class equationKCL : public equation {
public:
  equationKCL(int n);
  virtual void print();
   int mNode;
protected:
private:
};
#endif //EQUATIONKCL_H

