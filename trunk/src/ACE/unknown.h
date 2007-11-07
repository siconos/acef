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
  //mIndex in ----x'-------x--------Zs---------Zns--------
   int mIndex;
  //mDynIndex only for unknown in x: ----x'-------x--------Zs---------Zns--------
   int mDynIndex;
  //mIndexInVector index in x or Zs or Zns.
  int mIndexInVector;
  component * mComponent;
protected:
private:
};
#endif //UNKNOWN_H

