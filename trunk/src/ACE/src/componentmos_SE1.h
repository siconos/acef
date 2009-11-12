/**
MOS component:\n
Spice syntax: \n
.model mosP_Sah PMOS LEVEL=1 KP=3.24e-5 VT0=0.6\n
mosP 9 7 0 0 mosP_Sah \n


A model transistor is a necessary information. It specifies witch kind of transistor is used.\n
LEVEL=1 is a necessary option for the parser, but is not used in ACEF.\n
VT0 is the zero-bias voltage.\n
KP is the transconductance parameter.\n
The others parameters of the model are ignored.\n
ACEF use a piecewise linear model for the transistor. The number of hyperplan is a parameter available in the main program.\n

ace formulation:\n
\f[
dim(\lambda) = 2*NbHyperPlan
\f]
\f[
Zns = Ids
\f]
\f[
Zns = (mcoefs..-mcoefs)*\lambda
\f]
\f[
Y = I*\lambda + \left(\begin{array}{ccc}
0&-b&b\\
.&.&.\\
b&-b&0\\
.&.&.\\
\end{array}\right)*\left(\begin{array}{c}
Vs\\
Vg\\
Vd\\
\end{array}\right)
+Hyp
\f]


*/


/************************************************************************
  			componentmos_SE1.h 
**************************************************************************/

#ifndef COMPONENTMOS_SE1_H
#define COMPONENTMOS_SE1_H
#include "ace.h"

#include "component_linear_ns.h"

// Class componentMOS
// 
// 
class componentMOS_SE1 : public component_LINEAR_NS {
 dataMOS1 mData;
  int mNodeD;
  int mNodeG;
  int mNodeS;
  int mMode;
  int mNbHyp;
  ACE_DOUBLE mB;
  ACE_DOUBLE mVt;
  ACE_DOUBLE *mCoefs;
  ACE_DOUBLE *mHyp;
public:
  /*
    nbHyp is the number of hyperplan.
    if nbHyp==0, the model with 5 hp p 148 Bokhoven is used
  */
  componentMOS_SE1(dataMOS1 *d,int nbHyp);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  void stampMNA_V();
  virtual ~componentMOS_SE1();
  virtual void print();

protected:
private:
};
#endif //COMPONENTMOS_SE1_H

