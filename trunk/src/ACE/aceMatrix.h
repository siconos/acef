#ifndef ACEMATRIX_H
#define ACEMATRIX_H

#include "SimpleMatrix.h"

using namespace std;

class aceMatrix : public SimpleMatrix {
public:
  aceMatrix (unsigned int, unsigned int, UBLAS_TYPE=DENSE, unsigned int = 1, unsigned int = 1);

  virtual void display(ostream& os = cout) const;
  aceMatrix& operator = (const SimpleMatrix& );
  friend ostream & operator<<(ostream &f, const aceMatrix &Mat);


};

#endif //ACEMATRIX_H

