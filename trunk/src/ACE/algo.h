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
  algo(char * file);
  virtual ~algo();
  void perform();
  void printComponents();
  void stamp();
  void stampAfterInvertion();
  void preparStep();
  void simulate();

  components mInds;
  components mCaps;
  components mRess;
  components mIsrcs;
  components mVsrcs;
  components mDios;
  char mFile[ACE_CHAR_LENGTH];
  ofstream* mSimuStream;
  char mSimuFile[ACE_CHAR_LENGTH];

protected:
private:
    

};
#endif //ALGO_H

