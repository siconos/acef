#ifndef ACEMATRIX_H
#define ACEMATRIX_H

#include "SimpleMatrix.h"
#include "aceVector.h"
using namespace std;

class aceMatrix : public SimpleMatrix {
public:
  aceMatrix (unsigned int, unsigned int, UBLAS_TYPE=DENSE, unsigned int = 1, unsigned int = 1);

  virtual void display(ostream& os = cout) const;
  void setValueIfNotNull(unsigned int, unsigned int, double);
  MATRIX_UBLAS_TYPE getMat(){return mat;}

  aceMatrix& operator = (const SimpleMatrix& );
  friend ostream & operator<<(ostream &f, const aceMatrix &Mat);
  void MatrixToFortran(double * t);
  void FortranToMatrix(double * t);
  void MatrixToPath(int * I,int * J,double * t);
  void PathToMatrix(int * I,int * J,double * t);
  void set(const aceMatrix& M);
  static aceMatrix * load(char *file,UBLAS_TYPE typ);
  int getDimRow(){return dimRow;}
  int getDimCol(){return dimCol;}

  friend void ACEprod(const aceMatrix& A, const aceMatrix& B, aceMatrix& C, bool init = true);
  friend void ACEprod(const aceMatrix& A, const aceVector& x, aceVector& b, bool init = true);

};
TYPEDEF_SPTR(aceMatrix);

#endif //ACEMATRIX_H

