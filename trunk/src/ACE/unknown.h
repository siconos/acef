/************************************************************************
  			unknown.h 
**************************************************************************/

#ifndef UNKNOWN_H
#define UNKNOWN_H
#include "ace.h"

#include "component.h"
// Class unknown
// 
// 
class unknown {
public:
   unknown(int type, component *c);
   unknown(int type, int node);
  void print();
  void printdev();
  int mNode;
   int mType;
   int mIndex;
   int mDynIndex;
  component * mComponent;
protected:
private:
};
#endif //UNKNOWN_H

