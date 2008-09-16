


/** The component class.
*/


/**
 *
 *
 *Convention:\n
 *When possible, current go from Neg node to Pos node(for KCL law).\n
 *tension = V_neg - V_pos\n
 *
 * here a image\n
 * \image html componentConvention.jpg
**************************************************************************/

/*
 *Convetion:
 *
 *      UAB
 *   <---------- 
 *  A    ____    B
 *--|---|    |->-|---
 *       ----  I
 *UAB = VA-VB
 *
 *here:
 *       UNP
 *  <---------- 
 *  N    ____    P
 *--|---|    |->-|---
 *       ----  I
 *UNP = VN-VP
 *EXEMPLE
 *UNP = RI
 */
 
#ifndef COMPONENT_H
#define COMPONENT_H
#include "ace.h"
// Class component
// 
// 
#include "linearsystem.h"
class component {
public:
  component();
  virtual ~component();
  /**
   *\brief This methode fills the system:\n
   *Ax'=Bx+CZs+DZns+s\n
   * 0 =Ex+FZs+DZns+t
   */
  virtual void  stamp ();
  /**
   *\brief This methode fills the right side of the equations\n
   */
  virtual void stampTimer();
  /**
   *\brief Each component add its own unknowns.\n
   */
  virtual void  addUnknowns ();
  /**
   *\brief Each component add its constitutive equations.\n
   */
  virtual void  addEquations ();
  /**
   *\brief Display the component.\n
   */
  virtual void print ();

  /**
   * Null for the current defined branch.\n
   * Useful for the voltage defined component.\n
   *
   */
  unknown *mI;
  /**
   * An integer to defined the type of component.
   * For more details see \ref ace.h
   */
  int mType;
  /**
   * Id of the positive node.
   *
   */
  int mNodePos;
  /**
   * Id of the negative node.
   *
   */
  int mNodeNeg;
  unknown *mU;
  /**
   * Constitutive equation.
   *
   */
  equation *mEquation;
  /**
   * Name from the netlist.
   *
   */
  char *mName;
protected:
private:
    

};
#endif //COMPONENT_H

