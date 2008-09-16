/************************************************************************
  			componentcomp.h 

**************************************************************************/

#ifndef COMPONENTCOMP_H
#define COMPONENTCOMP_H
#include "ace.h"
#include "componentnlinear.h"

// Class componentCOMP
// 
// 
class componentCOMP : public componentNLINEAR {
public:
  dataCOMP mData;
  ACE_DOUBLE mV2;
  ACE_DOUBLE mV1;
  ACE_DOUBLE mEpsilon;
  ACE_DOUBLE mD11;
  ACE_DOUBLE mD12;
  int mNodeS;
  unknown *mVns;
  
  componentCOMP(dataCOMP *d);
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void  stamp ();
  virtual void  print ();

  virtual ~componentCOMP();
protected:
  
private:
 
};
#endif //COMPONENTCOMP_H

