


/** The  node.
*/


/**
 *
 *
**************************************************************************/

 
#ifndef NODE_H
#define NODE_H
#include "ace.h"
// Class node
// 
// 
class componentCAP;
#include "aceMatrix.h"
class node {
public:
  node(int i);
  virtual ~node();


  void stampV(ACE_DOUBLE v, ACE_DOUBLE * X, ACE_DOUBLE * Zs);
  void stampV(ACE_DOUBLE v,int lin,aceMatrix * X, aceMatrix * Zs);
  double getValue(aceVector * X, aceVector * Zs);
  void setIndexInZs(int i){mIndexInZs =i;};
  int getIndexInZs(){return mIndexInZs;};
  //  int getId(){return mId;}

  //  void alloc(int dimX,int dimZs);
  void setCapa(componentCAP *c){mC=c;};
  bool isVUnknown();
  void print();
  //  void VProd(ACE_DOUBLE * v);
  //  ACE_DOUBLE * getXBuffer(){return mXBuffer;};
  //  ACE_DOUBLE * getZsBuffer(){return mZsBuffer;};
protected:
//   int mDimX;
//   int mDimZs;
//   ACE_DOUBLE * mXCoef;
//   ACE_DOUBLE * mXBuffer;
//   ACE_DOUBLE * mZsCoef;
//   ACE_DOUBLE * mZsBuffer;
  componentCAP * mC;
  int mIndexInZs;
  int mId;
private:
    

};
#endif //NODE_H

