#ifndef ACEMATRIX_H
#define ACEMATRIX_H

#include "SimpleMatrix.h"

class aceMatrix : public SimpleMatrix {
public:
  aceMatrix (unsigned int, unsigned int, UBLAS_TYPE=DENSE, unsigned int = 1, unsigned int = 1);
  virtual void display();
  aceMatrix& operator = (const SimpleMatrix& );

};

#endif //ACEMATRIX_H

