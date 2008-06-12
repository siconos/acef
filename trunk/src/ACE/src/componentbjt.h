/************************************************************************
  			componentbjt.h

Bipolar Junction Transistor

      /c
 b  |/
----|
    |
    |\
      \e
       >

compementarity model:
---------------------

*)Amplification:

Ic=0+U-V
j=-alpha*(Vc-Ve) + U
k=V+betha*Ib-alpha*(Vc-Ve)

*)Condition Amplification:

Ib=W
l=-(Vb-Ve-Vbiais)

*)voltage repartition ???

Ub=Vb = Vc  + Vbiais -Z
m = Z - (Vc)

ace formulation:
----------------

lambda=(U,V,W)
Y=(j,k,l)


Zns = (Ic,Ib)

      |1 -1 0  |          |0|  | 0 0 0  |
Zns = |0 0  1  |*lambda + |0| +| 0 0 0  |(Vc,Vb,Ve)

   |1 0 0 |         |-alpha  0   alpha|           |0  0      |
Y= |0 1 0 |*lambda+ |-alpha  0   alpha|(Vc,Vb,Ve)+|0  betha  | Zns +CST
   |0 0 0 |         |0       -1  1    |           |0  0      |
   
CST=(0,0,Vbiais)

0<lambda per Y>0


**************************************************************************/

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

