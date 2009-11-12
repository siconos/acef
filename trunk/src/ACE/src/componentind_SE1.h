/************************************************************************
  			componentind_SE1.h 
**************************************************************************/

#ifndef COMPONENTIND_SE1_H
#define COMPONENTIND_SE1_H
#include "ace.h"
#include "componentdyn.h"

// Class componentIND_SE1
// 
// 
class componentIND_SE1 : public componentDYN {

public:
  dataIND mData;
  componentIND_SE1(dataIND * d);

  virtual void addUnknowns();
  virtual void addEquations();
  virtual void stamp();

virtual ~componentIND_SE1();

protected:
  
private:
  
};
#endif //COMPONENTIND_SE1_H

