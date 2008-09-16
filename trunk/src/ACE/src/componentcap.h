/************************************************************************
  			componentcap.h
**************************************************************************/

#ifndef COMPONENTCAP_H
#define COMPONENTCAP_H
#include "ace.h"

#include "componentdyn.h"
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
class componentCAP : public componentDYN {
public:
  dataCAP mData;
  componentCAP(dataCAP *d);
  ACE_DOUBLE * mICoefs;
  virtual ~componentCAP();
  /*
   *
   *
   */
  void  stampBeforeInvertion ();
  /*
   *
   *
   */
  void  stampAfterInvertion ();
  /*
   *
   *
   */
  virtual void addUnknowns();
  virtual void addEquations();

  /*
   *
   *
   */
  void addCurrentEquation();
  /*
   *
   *
   */
  void addTensionEquation();
  
  /**
   *Add the constitutive equation.
   *
   */
  void addCurrentUnknown();
  /**
   *Add Unp as an unknown.
   *
   */
  void addTensionUnknown();
  equationTEN* mTenEq;

  void printI();

protected:
 
private:
};
#endif //COMPONENTCAP_H



