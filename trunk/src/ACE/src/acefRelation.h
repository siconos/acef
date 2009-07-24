#ifndef ACEFRELATION_H
#define ACEFRELATION_H

#define WITH_KERNEL_RELATION

#ifdef WITH_KERNEL_RELATION

#include "SiconosKernel.hpp"
#include "ace.h"

class acefRelation : public FirstOrderType2R
{
protected:
  SimpleMatrix* mB;
  SimpleMatrix* mC;
  SimpleMatrix* mD;
  ACE_DOUBLE mLastComputedSource;
public:
  acefRelation();
  virtual ~acefRelation(){};


  virtual void initialize(SP::Interaction inter);
  void initJac(SimpleMatrix* C,SimpleMatrix* D,SimpleMatrix* B);

  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeH(double) ;
	
  /** default function to compute g 
   *  \param double : current time
   */
  virtual void computeG(double) ;

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJacH(double, unsigned int);
	
  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacG(double, unsigned int);
	


  double source(double t);

};

TYPEDEF_SPTR(acefRelation);
#endif

#endif
