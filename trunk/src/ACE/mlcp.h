/************************************************************************
  			mlcp.h 
**************************************************************************/

#ifndef MLCP_H
#define MLCP_H
#include "ace.h"
#include "aceMatrix.h"
#include <list>

#define SOLVER_ENUM 0
#define SOLVER_SIMPLEX 1

typedef std::list<unsigned long> ulongs;
typedef std::list<unsigned long>::iterator Itulongs;


// Class mlcp
// 
// 
class mlcp {
public:
  mlcp(unsigned int Dlcp,unsigned int Dlin);
  virtual ~mlcp();
  bool solve();

  void addGuess(unsigned long l);
  int mSolverType;
  
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
  bool mUseGuess;
  bool mTryGuess;
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
  bool solveWithSimplex();

};
#endif //MLCP_H
