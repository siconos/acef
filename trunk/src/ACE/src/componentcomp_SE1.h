/**

Comparator component.\n
ACEF syntax : .comp N+ N- Noutput Vmin=0 Vmax=3 Vepsilon=0.1\n

characteristic:\n


\image html comparator.jpg


acef formulation:\n
\f[
dim(\lambda) = 2
\f]
\f[
Zns=Vs
\f]
\f[
Y = \left(\begin{array}{c}
Y_1\\
Y_2\\
\end{array}\right)=

\left(\begin{array}{c}
Vp-Vn\\
Vp-Vn\\
\end{array}\right)+I\lambda+
\left(\begin{array}{c}
epsilon\\
-epsilon\\
\end{array}\right)
\f]
\f[
Zns = Vplus + (d11,d12)\lambda
\f]
*/

/************************************************************************
  			componentcomp_SE1.h 

**************************************************************************/

#ifndef COMPONENTCOMP_SE1_H
#define COMPONENTCOMP_SE1_H
#include "ace.h"
#include "component_linear_ns.h"

// Class componentCOMP
// 
// 
class componentCOMP_SE1 : public component_LINEAR_NS {
public:
  dataCOMP mData;
  ACE_DOUBLE mV2;
  ACE_DOUBLE mV1;
  ACE_DOUBLE mEpsilon;
  ACE_DOUBLE mD11;
  ACE_DOUBLE mD12;
  int mNodeS;
  unknown *mVns;
  
  componentCOMP_SE1(dataCOMP *d);
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void  stamp ();
  //  void stampMNA_V();

  virtual void  print ();

  virtual ~componentCOMP_SE1();
protected:
  
private:
 
};
#endif //COMPONENTCOMP_SE1_H

