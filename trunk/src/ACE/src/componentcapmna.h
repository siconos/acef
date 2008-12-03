/************************************************************************
  			componentcapmna.h
**************************************************************************/

#ifndef COMPONENTCAPMNA_H
#define COMPONENTCAPMNA_H
#include "ace.h"

#include "componentcap.h"
#include "equationten.h"

// Class componentCAP
// 
//
/**
 
 It implements the capacitor component.\n
 
 The constitutive equation is:\n
 \f[
 U_{np} = V_n-V_p
 \f]
 \f[
 C*\frac{U_{np}}{dt} = i(t)
 \f]
 
 */
class componentCAPMNA : public componentCAP {
public:
  componentCAPMNA(dataCAP *d);
  virtual ~componentCAPMNA();
  /*
   *
   *
   */
  virtual void  stampBeforeInvertion ();
  virtual void  stamp ();
  virtual void addUnknowns();
  virtual void addEquations();

protected:
 
private:
};
#endif //COMPONENTCAPMNA_H



