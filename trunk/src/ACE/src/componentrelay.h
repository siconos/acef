/************************************************************************
  			componentrelay.h 

**************************************************************************/

#ifndef COMPONENTRELAY_H
#define COMPONENTRELAY_H
#include "ace.h"
#include "componentnlinear.h"

// Class componentRELAY
// 
//
/**
   Electric relay:\n

   \image html componentrelay.jpg

   
   
   \image html componentrelayV1supV2.jpg
   ace formulation:\n
   \f[
 Y = \left(\begin{array}{c}
Z_1\\
W^-\\
\end{array}\right) \qquad
 L = \left(\begin{array}{c}
W^+\\
Z_2\\
\end{array}\right) 
   \f]
   \f[
   U_{ns} = V1-coef*Z_2
   \f]
   \f[
   Y=\left(\begin{array}{c}
Z_1\\
W^-\\
\end{array}\right) = \left(\begin{array}{c}
coef*(U_{ns} - V_2)\\
-U_e +W^++epsilon\\
\end{array}\right)
   \f]
\f[
0 \leq Y \, \perp \, lambda \geq 0
\f]

V1>V2 : coef=1\n
V2>V1 : coef=-1\n



 */
class componentRELAY : public componentNLINEAR {
public:
  dataCOMP mData;
  ACE_DOUBLE mV2;
  ACE_DOUBLE mV1;
  ACE_DOUBLE mOffset;
  int mNodeS;
  unknown *mVns;
  
  componentRELAY(dataCOMP *d);
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void  stamp ();
  virtual void  print ();

  virtual ~componentRELAY();
protected:
  
private:
 
};
#endif //COMPONENTRELAY_H

