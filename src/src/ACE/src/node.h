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

/*! \file node.h

*/



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

