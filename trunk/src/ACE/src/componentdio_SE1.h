/**

Diode component.\n
Spice syntax : D1 np nn DIODEF\n

acef formulation:


\f[
Zns=(Id)
\f]
\f[
Zns=\lambda
\f]

\f[
Y=Vn-Vp
\f]

\f[
0 \leq Y \, \perp \, \lambda \geq 0
\f]

*/
/************************************************************************
  			componentdio_SE1.h 
**************************************************************************/

#ifndef COMPONENTDIO_SE1_H
#define COMPONENTDIO_SE1_H
#include "ace.h"

#include "component_linear_ns.h"

// Class componentDIO
// 
// 
class componentDIO_SE1 : public component_LINEAR_NS {
 dataDIO mData;
  ACE_DOUBLE mThreshold;
public:
  componentDIO_SE1(dataDIO *d);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  virtual ~componentDIO_SE1();
  virtual void print ();
protected:
private:
};
#endif //COMPONENTDIO_SE1_H

