/************************************************************************
  			componentmos.h 
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
  int mN;
  ACE_DOUBLE mB;
  ACE_DOUBLE mVt;
  ACE_DOUBLE *mCoefs;
  ACE_DOUBLE *mHyp;
public:
  componentMOS(dataMOS1 *d);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  virtual ~componentMOS();
  virtual void print();

protected:
private:
};
#endif //COMPONENTMOS_H

