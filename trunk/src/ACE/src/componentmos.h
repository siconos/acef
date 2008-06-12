/************************************************************************
  			componentmos.h 

Zns = Ids
Zns = (mcoefs..-mcoefs)*lambda
            |0 -b b|
Y= I*lambda+|0 -b b|(Vs,Vg,Vd)+Hyp
            |b -b 0|
            |b -b 0|

mcoefs contains the variation of Ids gradient
Hyp contains the location of the variation of Ids gradient

**************************************************************************/

#ifndef COMPONENTMOS_H
#define COMPONENTMOS_H
#include "ace.h"

#include "componentnlinear.h"

// Class componentMOS
// 
// 
class componentMOS : public componentNLINEAR {
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
  componentMOS(dataMOS1 *d,int nbHyp);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  virtual ~componentMOS();
  virtual void print();

protected:
private:
};
#endif //COMPONENTMOS_H

