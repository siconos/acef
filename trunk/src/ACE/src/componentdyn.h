/************************************************************************
  			componentdyn.h 
**************************************************************************/

#ifndef COMPONENTDYN_H
#define COMPONENTDYN_H
#include "ace.h"

#include "componentlinear.h"

/** The componentDYN class.
*/


/**
 *
 *It has a specufic equations menber.\n
 *
**************************************************************************/
class componentDYN : public componentLINEAR {
public:
  componentDYN();
  virtual ~componentDYN();
  equation *mDynEquation;
protected:
private:
};
#endif //COMPONENTDYN_H

