/*

Bipolar Junction Transistor component.\n
Spice syntax : Q23 nc nb ne ...\n
*/
/*      /c
   b  |/
  ----|
      |
      |\
        \e
         >
*/
/**

compementarity model:\n
---------------------\n

*)Amplification:\n
\f[
I_c=0+U-V
\f]
\f[
j=-alpha*(Vc-Ve) + U
\f]
\f[
k=V+betha*Ib-alpha*(Vc-Ve)
\f]
*)Condition Amplification:\n
\f[
Ib=W
\f]

\f[
l=-(Vb-Ve-Vbiais)
\f]

*)voltage repartition ???
\f[

Ub=Vb = Vc  + Vbiais -Z
\f]
\f[
m = Z - (Vc)
\f]

ace formulation:
----------------
\f[
lambda=(U,V,W)
\f]
\f[
Y=(j,k,l)
\f]

\f[
Zns = (Ic,Ib)
\f]

\f[
Zns = \left(\begin{array}{ccc}
1&-1&0\\
0&0&1\\
\end{array}\right)
*lambda +
\left(\begin{array}{c}
0\\
0\\
\end{array}\right)
+
\left(\begin{array}{ccc}
0&0&0\\
0&0&0\\
\end{array}\right)*
\left(\begin{array}{c}
V_c\\
V_b\\
V_e\\
\end{array}\right)

\f]
\f[
Y = \left(\begin{array}{ccc}
1&0&0\\
0&0&1\\
0&0&0\\
\end{array}\right)
*lambda +
\left(\begin{array}{ccc}
-alpha&0&alpha\\
-alpha&0&alpha\\
0&-1&1\\
\end{array}\right)*
\left(\begin{array}{c}
V_c\\
V_b\\
V_e\\
\end{array}\right)
+
\left(\begin{array}{cc}
0&0\\
0&betha\\
0&0\\
\end{array}\right)*Zns
+CST
\f]

\f[
CST=(0,0,Vbiais)
\f]
\f[
0 \leq Y \, \perp \, lambda \geq 0
\f]


**************************************************************************/
/*
      |1 -1 0  |          |0|  | 0 0 0  |
Zns = |0 0  1  |*lambda + |0| +| 0 0 0  |(Vc,Vb,Ve)

   |1 0 0 |         |-alpha  0   alpha|           |0  0      |
Y= |0 1 0 |*lambda+ |-alpha  0   alpha|(Vc,Vb,Ve)+|0  betha  | Zns +CST
   |0 0 0 |         |0       -1  1    |           |0  0      |
*/
#ifndef COMPONENTBJT_H
#define COMPONENTBJT_H
#include "ace.h"

#include "componentnlinear.h"

// Class componentBJT
// 
// 
class componentBJT : public componentNLINEAR {
 dataBJT mData;
  int mNodeB;
  int mNodeC;
  int mNodeE;
  int mMode;
  ACE_DOUBLE mBetha;/*current ampli*/
  ACE_DOUBLE mVBiais;/*minimal value for amplification*/
  ACE_DOUBLE mAlpha;
public:
  /*
  */
  unknown *mIb;
  unknown *mIc;

  componentBJT(dataBJT *d);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  virtual ~componentBJT();
  virtual void print();

protected:
private:
};
#endif //COMPONENTBJT_H

