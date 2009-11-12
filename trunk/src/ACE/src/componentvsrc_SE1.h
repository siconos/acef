/************************************************************************
  			componentvsrc_SE1.h 

**************************************************************************/

#ifndef COMPONENTVSRC_SE1_H
#define COMPONENTVSRC_SE1_H
#include "componentlinear.h"
#include "ace.h"
// Class componentVSRC
// 
// 
class componentVSRC_SE1 : public componentLINEAR {

public:
  dataVSRC mData;
  ACE_DOUBLE mCurrentValue;
  componentVSRC_SE1(dataVSRC *d);
  virtual ~componentVSRC_SE1();
  virtual void  stamp ();
  virtual void stampTimer();
  void stampTime();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  double * mpCurValue;
protected:
private:
};
#endif //COMPONENTVSRC_SE1_H

