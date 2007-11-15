/************************************************************************
  			mlcp.h 
**************************************************************************/

#ifndef MLCP_H
#define MLCP_H
#include "ace.h"
#include "aceMatrix.h"

// Class mlcp
// 
// 
class mlcp {
public:
  mlcp(unsigned int Dlcp,unsigned int Dlin);
  virtual ~mlcp();
  bool solve();
  
  aceMatrix *mW1;
  aceMatrix *mZ1;
  aceMatrix *mZ2;
  int* mW1Z1;
  unsigned int mDlcp;
  unsigned int mDlin;
  unsigned long mCurEnum;
  unsigned long mCmp;
  unsigned long mMaxEnum;
  
  
  aceMatrix *mQ1;
  aceMatrix *mQ2;
  aceMatrix *mQ;
  aceMatrix *mM11;
  aceMatrix *mM12;
  aceMatrix *mM21;
  aceMatrix *mM22;
  aceMatrix *mM;

  void printInPut(ostream& os = cout);
  void printOutPut(ostream& os = cout);
protected:
private:
  void initEnum();
  bool nextEnum();
};
#endif //MLCP_H
