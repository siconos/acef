/************************************************************************
  			algo.h 

**************************************************************************/
/**The algo class.
 */

/**
*/
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
  static algo *sAlgo;
  /*
   * file : a netlist.
   */
  algo(char * file);
  virtual ~algo();
  void perform();
  void printComponents();
  void stamp();
  void stampAfterInvertion();
  void preparStep(double time);
  void simulate();
  void parseComponents();
  void parseComponents_SE();
  
  void computeNonLinearEquations(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosVector& SICONOS_H); 
  void computeNonLinearJacL_H(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosMatrix& SICONOS_D);




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
  component_NLs mMos_NL;
  components mBjt;
  components mComps;
  components mRelays;
  char mFile[ACE_CHAR_LENGTH];
  ofstream* mSimuStream;
  char mSimuFile[ACE_CHAR_LENGTH];

protected:
  void performSemiExplicit();
  void performWithOutInvert();
  void performMNA();
  void performMNA_V();
  void performSE();
private:
    

};
#endif //ALGO_H

