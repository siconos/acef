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
  void MatrixToFortran(double * t);
  void FortranToMatrix(double * t);
  void MatrixToPath(int * I,int * J,double * t);
  void PathToMatrix(int * I,int * J,double * t);


};

#endif //ACEMATRIX_H

