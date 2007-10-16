/************************************************************************
  			componentisrc.h 

**************************************************************************/

#ifndef COMPONENTISRC_H
#define COMPONENTISRC_H
#include "ace.h"
#include "componentlinear.h"

// Class componentISRC
// 
// 
class componentISRC : public componentLINEAR {

public:
  dataISRC mData;
  componentISRC(dataISRC *d);
  virtual void stamp();
  virtual void print();
protected:
private:
};
#endif //COMPONENTISRC_H

