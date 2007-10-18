/************************************************************************
  			componentcap.h
**************************************************************************/

#ifndef COMPONENTCAP_H
#define COMPONENTCAP_H
#include "ace.h"

#include "componentdyn.h"
#include "equationten.h"

// Class componentCAP
// 
// 
class componentCAP : public componentDYN {
public:
  dataCAP mData;
  componentCAP(dataCAP *d);
  ACE_DOUBLE * mICoefs;
  virtual ~componentCAP();
  void  stampBeforeInvertion ();
  void  stampAfterInvertion ();
  virtual void addUnknowns();
  virtual void addEquations();

  void addCurrentEquation();
  void addTensionEquation();
  
  void addCurrentUnknown();
  void addTensionUnknown();
  equationTEN* mTenEq;
protected:
 
private:
};
#endif //COMPONENTCAP_H



