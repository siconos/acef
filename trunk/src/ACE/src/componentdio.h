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
  ACE_DOUBLE mThreshold;
public:
  componentDIO(dataDIO *d);
  virtual void  addUnknowns ();
  virtual void  stamp ();
  virtual ~componentDIO();
  virtual void print ();
protected:
private:
};
#endif //COMPONENTDIO_H

