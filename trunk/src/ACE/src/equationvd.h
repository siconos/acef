/************************************************************************
  			equationvd.h 
**************************************************************************/

#ifndef EQUATIONVD_H
#define EQUATIONVD_H
#include "equation.h"
#include "ace.h"

class equationVD : public equation {
public:
  equationVD(char* name);
  virtual void print();
protected:
private:
};
#endif //EQUATIONVD_H

