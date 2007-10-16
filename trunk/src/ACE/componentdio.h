/************************************************************************
  			componentdio.h 
**************************************************************************/

#ifndef COMPONENTDIO_H
#define COMPONENTDIO_H
#include "ace.h"

#include "componentnlinear.h"

// Class componentDIO
// 
// 
class componentDIO : public componentNLINEAR {
 dataDIO mData;
public:
  componentDIO(dataDIO *d);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  virtual ~componentDIO();
protected:
private:
};
#endif //COMPONENTDIO_H

