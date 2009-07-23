/************************************************************************
  			component_none_linear_ns.h 
**************************************************************************/

#ifndef COMPONENT_NONE_LINEAR_NS_H
#define COMPONENT_NONE_LINEAR_NS_H
#include "ace.h"
#include "component_ns.h"

// Class component_NONE_LINEAR_NS
// 
// 

class component_NONE_LINEAR_NS : public component_NS {
  public:
  component_NONE_LINEAR_NS():component_NS()
  {}
  virtual void computeNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosVector& SICONOS_H){;}
  virtual void computeJacNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosMatrix& SICONOS_D){;}
};

#endif //COMPONENT_NONE_LINEAR_NS_H

