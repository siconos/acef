/************************************************************************
  			algo.h 

**************************************************************************/

#ifndef ALGO_H
#define ALGO_H
#include "ace.h"
#include "linearsystem.h"

// Class algo
// 
// 
class algo {
public:
  static linearSystem sls;
  algo();
  virtual ~algo();
  void perform();
  void printComponents();
  void stamp();

  components mInds;
  components mCaps;
  components mRess;
  components mIsrcs;
  components mVsrcs;
  components mDios;
  
protected:
private:
    

};
#endif //ALGO_H

