/************************************************************************
  			componentvcvs.h 

**************************************************************************/

#ifndef COMPONENTVCVS_H
#define COMPONENTVCVS_H
#include "componentlinear.h"
#include "ace.h"
// Class componentVCVS
// 
// 
class componentVCVS : public componentLINEAR {

public:
  dataVCVS mData;
  ACE_DOUBLE mCurrentValue;
  componentVCVS(dataVCVS *d);
  virtual ~componentVCVS();
  virtual void  stamp ();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void print();
protected:
private:
};
#endif //COMPONENTVCVS_H

