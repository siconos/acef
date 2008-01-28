#ifndef ACEVECTOR_H
#define ACEVECTOR_H

#include "SimpleVector.h"

using namespace std;

class aceVector : public SimpleVector {
public:
  aceVector (unsigned int,  UBLAS_TYPE=DENSE);

  virtual void display(ostream& os = cout) const;
  void setValueIfNotNull(unsigned int, double);

  aceVector& operator = (const SimpleVector& );
  friend ostream & operator<<(ostream &f, const aceVector &Mat);
  void VectorToFortran(double * t);
  void FortranToVector(double * t);
  void VectorToPath(int * I,double * t);
  void PathToVector(int * I,double * t);
  void set(const aceVector& V);
  unsigned int dimRow;

};

#endif //ACEVECTOR_H

