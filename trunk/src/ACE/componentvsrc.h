/************************************************************************
  			componentvsrc.h 

**************************************************************************/

#ifndef COMPONENTVSRC_H
#define COMPONENTVSRC_H
#include "componentlinear.h"
#include "ace.h"
// Class componentVSRC
// 
// 
class componentVSRC : public componentLINEAR {

public:
  dataVSRC mData;
  ACE_DOUBLE mCurrentValue;
  componentVSRC(dataVSRC *d);
  virtual ~componentVSRC();
  virtual void  stamp ();
  void stampTime();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
protected:
private:
};
#endif //COMPONENTVSRC_H

