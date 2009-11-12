/************************************************************************
  			componentvcvs_SE1.h 

**************************************************************************/

#ifndef COMPONENTVCVS_SE1_H
#define COMPONENTVCVS_SE1_H
#include "componentlinear.h"
#include "ace.h"
// Class componentVCVS_SE1
// 
// 
class componentVCVS_SE1 : public componentLINEAR {

public:
  dataVCVS mData;
  ACE_DOUBLE mCurrentValue;
  componentVCVS_SE1(dataVCVS *d);
  virtual ~componentVCVS_SE1();
  virtual void  stamp ();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void print();
protected:
private:
};
#endif //COMPONENTVCVS_SE1_H

