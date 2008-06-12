/************************************************************************
  			componentvccs.h 

**************************************************************************/

#ifndef COMPONENTVCCS_H
#define COMPONENTVCCS_H
#include "componentlinear.h"
#include "ace.h"
// Class componentVCCS
// 
// 
class componentVCCS : public componentLINEAR {

public:
  dataVCCS mData;
  componentVCCS(dataVCCS *d);
  virtual ~componentVCCS();
  virtual void  stamp ();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void print();
protected:
private:
};
#endif //COMPONENTVCCS_H

