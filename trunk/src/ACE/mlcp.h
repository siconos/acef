/************************************************************************
  			mlcp.h 
**************************************************************************/

#ifndef MLCP_H
#define MLCP_H
#include "ace.h"
#include "aceMatrix.h"
#include <list>


typedef std::list<unsigned long> ulongs;
typedef std::list<unsigned long>::iterator Itulongs;


// Class mlcp
// 
// 
class mlcp {
public:
  mlcp(unsigned int Dlcp,unsigned int Dlin,int solverType = ACE_SOLVER_ENUM);
  virtual ~mlcp();
  bool solve();
  void stopSolver();
  bool initSolver();
  void addGuess(aceMatrix *mZ);

  void addGuess(unsigned long l);
  void setCurrentConfig(unsigned long l);
  aceMatrix *mW1;
  aceMatrix *mZ1;
  aceMatrix *mZ2;
  int* mW1Z1;
  unsigned int mDlcp;
  unsigned int mDlin;
  unsigned long mCase;
  unsigned long mCurEnum;
  ulongs mGuess;
  Itulongs mItGuess;
  
  //CONFIGURATION--OPTIONS
  bool mUseGuess;
  int mSolverType;

  unsigned long mCmp;
  unsigned long mMaxEnum;
  double mPourCent;
  
  
  aceMatrix *mQ1;
  aceMatrix *mQ2;
  aceMatrix *mQ;
  aceMatrix *mM11;
  aceMatrix *mM12;
  aceMatrix *mM21;
  aceMatrix *mM22;
  aceMatrix *mM;

  void printGuess(ostream& os = cout);
  void printInPut(ostream& os = cout);
  void printInPutABCDab(ostream& os = cout);

  void printOutPut(ostream& os = cout);
protected:
private:
  void initEnum();
  bool nextEnum();
  bool tryGuess();
  void initGuess();
  void affectW1Z1(unsigned long ll);
  bool solveGuessAndIt();
  bool solveWithSimplex();
  bool solveWithPath();
  //INTERNAL--OPTION
  bool mTringGuess;
  bool mTryOnlyGuess;

  //
  double* mA;
  double* mB;
  double* mC;
  double* mD;
  double* ma;
  double* mb;
  double* mu;
  double* mv;
  double* mw;


};
#endif //MLCP_H
