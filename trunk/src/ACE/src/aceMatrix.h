/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file aceMatrix.h

*/
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
  //  void MatrixToPath(int * I,int * J,double * t);
  //  void PathToMatrix(int * I,int * J,double * t);
  void set(const aceMatrix& M);
  static aceMatrix * load(char *file,UBLAS_TYPE typ);
  int getDimRow(){return dimRow;}
  int getDimCol(){return dimCol;}
  void incValue(int lin, int col, double val);

  //  static void siconosProd(const SiconosMatrix& a, const SiconosVector& v1, SiconosVector& v2, bool b= true);

  //  friend void ACEprod(const aceMatrix& A, const aceMatrix& B, aceMatrix& C, bool init = true);
  //  friend void ACEprod(const aceMatrix& A, const aceVector& x, aceVector& b, bool init = true);
  friend void ACEprod(const SimpleMatrix& A, const SimpleMatrix& B, SimpleMatrix& C, bool init = true);
  friend void ACEprod(const SimpleMatrix& A, const SimpleVector& x, SimpleVector& b, bool init = true);
  friend void ACEbidonF();

};
TYPEDEF_SPTR(aceMatrix);

#endif //ACEMATRIX_H

