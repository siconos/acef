/************************************************************************
  			algo.h 

**************************************************************************/

#ifndef ALGO_H
#define ALGO_H
#include "ace.h"

// Class algo
// 
//
class linearSystem;
class algo {
public:
  static linearSystem *spls;
  algo(char * file);
  virtual ~algo();
  void perform();
  void printComponents();
  void stamp();
  void stampAfterInvertion();
  void preparStep(double time);
  void simulate();

  components mInds;
  components mCaps;
  components mRess;
  components mIsrcs;
  components mVsrcs;
  components mVcvs;
  components mVccs;
  components mArbs;
  components mDios;
  components mMos;
  components mBjt;
  components mComps;
  char mFile[ACE_CHAR_LENGTH];
  ofstream* mSimuStream;
  char mSimuFile[ACE_CHAR_LENGTH];

protected:
private:
    

};
#endif //ALGO_H

