/************************************************************************
  			componentcap.h
**************************************************************************/

#ifndef COMPONENTCAP_H
#define COMPONENTCAP_H
#include "ace.h"

#include "componentdyn.h"

// Class componentCAP
// 
// 
class componentCAP : public componentDYN {
public:
  dataCAP mData;
  componentCAP(dataCAP *d);
  virtual ~componentCAP();
  void  stampBeforeInvertion ();
  void  stampAfterInvertion ();
  virtual void addUnknowns();
  virtual void addEquations();
  
  void addCurrentUnknown();

protected:
 
private:
};
#endif //COMPONENTCAP_H



