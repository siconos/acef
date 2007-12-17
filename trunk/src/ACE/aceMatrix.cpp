#include "aceMatrix.h"



aceMatrix::aceMatrix ( unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
SimpleMatrix(row,col,typ,upper,lower)
{
  ;
}


void aceMatrix::display(ostream& os) const{
  os <<"["<<dimRow<<","<<dimCol<<"]"<<endl;
  for (unsigned int i=0;i<dimRow;i++){
    for (unsigned int j=0;j<dimCol;j++){
      os <<"\t"<<((SimpleMatrix&)*this)(i,j);
    }
    os <<"\n";
  }
      
}

aceMatrix& aceMatrix::operator = (const SimpleMatrix& m)
{
  ((SimpleMatrix&)*this) =  m;
  return *this;
}


ostream & operator<<(ostream &f, const aceMatrix &Mat)
{
  Mat.display(f);
  return f;
}

void aceMatrix::MatrixToFortran(double * t){
  for (unsigned int i=0; i < dimRow; i++)
    for (unsigned int j=0; j < dimCol; j++)
      t[j*dimRow+i]=((SimpleMatrix&)*this)(i,j);
    
}
void aceMatrix::FortranToMatrix(double * t){
  for (unsigned int i=0; i < dimRow; i++)
    for (unsigned int j=0; j < dimCol; j++){
      setValue(i,j,t[j*dimRow+i]);
    }
}
