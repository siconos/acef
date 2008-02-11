#include "aceMatrix.h"
#include "ace.h"
#include <fstream>


aceMatrix::aceMatrix ( unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
SimpleMatrix(row,col,typ,upper,lower)
{
  ;
}

void aceMatrix::set(const aceMatrix& M){
  unsigned int lin,col;
  ACE_times[ACE_TIMER_TEST].start();
  zero();
  for (lin = 0; lin < dimRow;lin++){
    for (col = 0; col < dimCol;col++){
      double v = M.getValue(lin,col);
      if (fabs(v) > ACE_NULL_COEF_MAT){
	setValue(lin,col,v);
      }
    }
  }
  ACE_times[ACE_TIMER_TEST].stop();
}
void aceMatrix::setValueIfNotNull(unsigned int lin, unsigned int col , double v){
  if (fabs(v) > ACE_NULL_COEF_MAT){
    setValue(lin,col,v);
  }else{
    if (num==1)
      setValue(lin,col,0);
    else{
      double *d=(*mat.Sparse).find_element(lin,col);
      if (d)
	(*mat.Sparse).erase_element(lin,col);
    }
  }
}

void aceMatrix::display(ostream& os) const{
  if (num == 1){
    os <<"DENSE MATRIX ["<<dimRow<<","<<dimCol<<"]"<<endl;
    for (unsigned int i=0;i<dimRow;i++){
      for (unsigned int j=0;j<dimCol;j++){
	os <<"\t"<<getValue(i,j);
      }
      os <<"\n";
    }
  }else if (num ==4){
    os <<"SPARSE MATRIX ["<<dimRow<<","<<dimCol<<"]"<<endl;
    for (unsigned int i=0;i<dimRow;i++){
      for (unsigned int j=0;j<dimCol;j++){
	double *d=(*mat.Sparse).find_element(i,j);
	if (d)
	  os <<"\t"<<*d;
	else
	  os <<"\t"<<"N";
      }
      os <<"\n";
    }
  }else
    os<<"unknown matrix type !!!!!\n";
      
}
aceMatrix * aceMatrix::load(char * file,UBLAS_TYPE typ){
  int i,j;
  aceMatrix *res=0;
  double aux;
  try{
    ifstream pin(file);
    ACE_CHECK_IERROR(pin,"aceMatrix::load error no file");
    pin >>i;
    pin >>j;
    res = new aceMatrix(i,j,typ);
    for (int ii=0;ii<i;ii++)
      for (int jj=0;jj<j;jj++){
	pin>>aux;
	res->setValueIfNotNull(ii,jj,aux);
      }
  }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("aceMatrix::load");
    }
  return res;
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
      t[j*dimRow+i]=getValue(i,j);
    
}
void aceMatrix::FortranToMatrix(double * t){
  for (unsigned int i=0; i < dimRow; i++)
    for (unsigned int j=0; j < dimCol; j++){
      setValue(i,j,t[j*dimRow+i]);
    }
}
void aceMatrix::MatrixToPath(int * I,int * J,double * t){
  ;
}
void aceMatrix::PathToMatrix(int * I,int * J,double * t){
  ;
}
void ACEprod(const aceMatrix& A, const aceMatrix& B, aceMatrix& C, bool init)
{
//   if (init && ACE_MAT_TYPE == DENSE && &A!=&C && &B!=&C)
//     gemm(A,B,C);
//   else
    prod(A,B,C,init);
}
void ACEprod(const aceMatrix& A, const aceVector& x, aceVector& b, bool init)
{
//   if (init && ACE_MAT_TYPE == DENSE)
//     gemv(A,x,b);
//   else
    prod(A,x,b,init);
}
