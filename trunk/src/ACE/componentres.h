/************************************************************************
  			componentres.h 

**************************************************************************/

#ifndef COMPONENTRES_H
#define COMPONENTRES_H
#include "ace.h"

#include "componentlinear.h"



// Class componentRES
// 
// 
class componentRES : public componentLINEAR {

public:
  dataRES mData;
  componentRES(dataRES *d);
  virtual void stamp();
  virtual void print();
  virtual ~componentRES();

protected:
private:
};
#endif //COMPONENTRES_H

