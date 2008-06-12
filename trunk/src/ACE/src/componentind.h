/************************************************************************
  			componentind.h 
**************************************************************************/

#ifndef COMPONENTIND_H
#define COMPONENTIND_H
#include "ace.h"
#include "componentdyn.h"

// Class componentIND
// 
// 
class componentIND : public componentDYN {

public:
  dataIND mData;
  componentIND(dataIND * d);

  virtual void addUnknowns();
  virtual void addEquations();
  virtual void stamp();

virtual ~componentIND();

protected:
  
private:
  
};
#endif //COMPONENTIND_H

