#include "aceVector.h"
#include "ace.h"


aceVector::aceVector ( unsigned int row, UBLAS_TYPE typ):
SimpleVector(row,typ)
{
 dimRow = row ;
}

void aceVector::set(const aceVector& V){
  unsigned int lin,col;
  ACE_times[ACE_TIMER_TEST].start();
  zero();
  for (lin = 0; lin < dimRow;lin++){
    double v = V.getValue(lin);
    if (fabs(v) > ACE_NULL_COEF_MAT){
      setValue(lin,v);
    }
  }
  ACE_times[ACE_TIMER_TEST].stop();
}
void aceVector::setValueIfNotNull(unsigned int lin,  double v){
  if (fabs(v) > ACE_NULL_COEF_MAT){
    setValue(lin,v);
  }else{
    if (num==1)
      setValue(lin,0);
    else{
      double *d=(*vect.Sparse).find_element(lin);
      if (d)
	(*vect.Sparse).erase_element(lin);
    }
  }
}

void aceVector::display(ostream& os) const{
  if (num == 1){
    os <<"DENSE VECTOR ["<<dimRow<<"]"<<endl;
    for (unsigned int i=0;i<dimRow;i++){
	os <<"\t"<<getValue(i);
	os <<"\n";
      }
  }else if (num ==4){
    os <<"SPARSE VECTOR ["<<dimRow<<"]"<<endl;
    for (unsigned int i=0;i<dimRow;i++){
      double *d=(*vect.Sparse).find_element(i);
      if (d)
	os <<"\t"<<*d;
      else
	os <<"\t"<<"N";
      os <<"\n";
    }
  }else
    os<<"unknown vector type !!!!!\n";
      
}

aceVector& aceVector::operator = (const SimpleVector& m)
{
  ((SimpleVector&)*this) =  m;
  return *this;
}


ostream & operator<<(ostream &f, const aceVector &Vec)
{
  Vec.display(f);
  return f;
}

void aceVector::VectorToFortran(double * t){
  for (unsigned int i=0; i < dimRow; i++)
      t[i]=getValue(i);
    
}
void aceVector::FortranToVector(double * t){
  for (unsigned int i=0; i < dimRow; i++)
      setValue(i,t[i]);
}
void aceVector::VectorToPath(int * I,double * t){
  ;
}
void aceVector::PathToVector(int * I,double * t){
  ;
}
