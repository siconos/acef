/************************************************************************
  			componentarb.h 

**************************************************************************/

#ifndef COMPONENTARB_H
#define COMPONENTARB_H
#include "componentlinear.h"
#include "ace.h"
// Class componentARB
// 
// 
class componentARB : public componentLINEAR {

public:
  dataARB mData;
  ACE_DOUBLE mCurrentValue;
  componentARB(dataARB *d);
  virtual ~componentARB();
  virtual void  stamp ();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
protected:
private:
};
#endif //COMPONENTARB_H

