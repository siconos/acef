/************************************************************************
  			componentdyn.h 
**************************************************************************/

#ifndef COMPONENTDYN_H
#define COMPONENTDYN_H
#include "ace.h"

#include "componentlinear.h"

// Class componentDYN
// 
// 
class componentDYN : public componentLINEAR {
public:
  componentDYN();
  virtual ~componentDYN();
  equation *mDynEquation;
protected:
private:
};
#endif //COMPONENTDYN_H

