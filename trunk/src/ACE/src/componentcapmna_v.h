/************************************************************************
  			componentcapmna_v.h
**************************************************************************/

#ifndef COMPONENTCAPMNA_V_H
#define COMPONENTCAPMNA_V_H
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
class componentCAPMNA_V : public componentCAP {
public:
  componentCAPMNA_V(dataCAP *d);
  virtual ~componentCAPMNA_V();
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
#endif //COMPONENTCAPMNA_V_H



