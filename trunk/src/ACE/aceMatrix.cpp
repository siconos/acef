#include "aceMatrix.h"



aceMatrix::aceMatrix ( unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
SimpleMatrix(row,col,typ,upper,lower)
{
  ;
}


void aceMatrix::display(){
  printf("[%d,%d]\n",dimRow,dimCol);
  for (int i=0;i<dimRow;i++){
    for (int j=0;j<dimCol;j++){
      printf("\t%f",((SimpleMatrix&)*this)(i,j));
    }
    printf("\n");
  }
      
}

aceMatrix& aceMatrix::operator = (const SimpleMatrix& m)
{
  ((SimpleMatrix&)*this) =  m;
  return *this;
}
