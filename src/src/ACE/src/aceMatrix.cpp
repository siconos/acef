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

/*! \file aceMatrix.cpp

*/
#include "aceMatrix.h"
#include "ace.h"
#include <fstream>

// void aceMatrix::siconosProd(const SiconosMatrix& a, const SiconosVector& v1, SiconosVector& v2, bool b)
// {
//   //prod(a,v1,v2,b);
// }
aceMatrix::aceMatrix ( unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
SimpleMatrix(row,col,typ,upper,lower)
{
  ;
}
void aceMatrix::incValue(int lin, int col, double val){
  setValue(lin,col,getValue(lin,col)+val);
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
	//	os <<"\t"<<getValue(i,j);
	if (j)
	  printf(",");
	printf(" %.10e",getValue(i,j));
      }
      os <<";\n";
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
// void aceMatrix::MatrixToPath(int * I,int * J,double * t){
//   ;
// }
// void aceMatrix::PathToMatrix(int * I,int * J,double * t){
//   ;
// }
void ACEprod(const SimpleMatrix& A, const SimpleMatrix& B, SimpleMatrix& C, bool init)
{
//   if (init && ACE_MAT_TYPE == DENSE && &A!=&C && &B!=&C)
//     gemm(A,B,C);
//   else

  //ACE_times[ACE_TIMER_PROD_MAT].start();

    prod(A,B,C,init);
    //ACE_times[ACE_TIMER_PROD_MAT].stop();
}
void ACEprod(const SimpleMatrix& A, const SimpleVector& x, SimpleVector& b, bool init)
{
//   if (init && ACE_MAT_TYPE == DENSE)
//     gemv(A,x,b);
//   else


  //ACE_times[ACE_TIMER_PROD_MAT].start();
    prod(A,x,b,init);
    //    ACE_times[ACE_TIMER_PROD_MAT].stop();

}
void ACEbidonF(){
  printf("a");
}
