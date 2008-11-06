/************************************************************************
  			mlcp.h
|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0


{M11 M12} Z1  a  W1
{       }*  +  =
{M21 M22} Z2  b  0

**************************************************************************/

#ifndef MLCP_ACE_H
#define MLCP_ACE_H
#include "ace.h"
#include "aceMatrix.h"
#include "aceVector.h"
#include <list>
#include "NonSmoothDrivers.h"


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
  void update();
  //  void addGuess(aceVector *mZ);

  //  void addGuess(unsigned long l);
  //  void setCurrentConfig(unsigned long l);
  aceVector *mW1;
  aceVector *mZ1;
  aceVector *mZ2;
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
  
  
  aceVector *mQ1;
  aceVector *mQ2;
  aceMatrix *mM11;
  aceMatrix *mM12;
  aceMatrix *mM21;
  aceMatrix *mM22;
  aceMatrix *mM;
  bool mTryM;
  
  //  void printGuess(ostream& os = cout);
  void printInPut(ostream& os = cout);
  void printInPutABCDab(ostream& os = cout);

  void printOutPut(ostream& os = cout);


  //Numerics interface
  MixedLinearComplementarity_Problem mProblem;
  Solver_Options mOptions;
  Numerics_Options mNumericsOptions;

protected:
private:
  //  void   fillSolution();

//   void initEnum();
//   bool nextEnum();
//   bool tryGuess();
//   void initGuess();
//   void affectW1Z1(unsigned long ll);
//   bool solveGuessAndIt();
//   bool solveWithSimplex();
//   bool solveWithPath();
  bool solveWithNumerics();
  bool solveLinearSystem();

  //INTERNAL--OPTION
  bool mTringGuess;
  bool mTryOnlyGuess;

  //
  double* mMd;
  //  double* mA;
  //  double* mB;
  //  double* mC;
  //  double* mD;
  //  double* ma;
  //  double* mb;
    double* mu;
    double* mv;
    double* mw;

  Index dim;
  Index start;
  Index end;

  int * mIwork;
  double * mDwork;


};
#endif //MLCP_ACE_H
