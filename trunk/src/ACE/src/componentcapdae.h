/************************************************************************
  			componentcapdae.h
**************************************************************************/

#ifndef COMPONENTCAPDAE_H
#define COMPONENTCAPDAE_H
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
class componentCAPDAE : public componentCAP {
public:
  componentCAPDAE(dataCAP *d);
  virtual ~componentCAPDAE();
  /*
   *
   *
   */
  virtual void  stampBeforeInvertion ();
  virtual void  stamp ();

protected:
 
private:
};
#endif //COMPONENTCAPDAE_H



