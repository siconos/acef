/************************************************************************
  			component.h 

Convention:
When possible, current go from Neg node to Pos node(for KCL law).
tension = V_neg - V_pos
stamp() to write : Ax'=Bx+CZs+DZns+s
                    0 =Ex+FZs+DZns+t
**************************************************************************/

#ifndef COMPONENT_H
#define COMPONENT_H
#include "ace.h"
// Class component
// 
// 
class component {
public:
  component();
  virtual ~component();
  virtual void  stamp ();
  virtual void stampTimer();
  virtual void  addUnknowns ();
  virtual void  addEquations ();
  virtual void print ();

  unknown *mI;
  int mType;
  int mNodePos;
  int mNodeNeg;
  unknown *mU;
  equation *mEquation;
  char *mName;
protected:
private:
    

};
#endif //COMPONENT_H

