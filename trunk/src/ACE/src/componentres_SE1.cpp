/************************************************************************
  			componentres_SE1.cpp 
**************************************************************************/

#include "componentres_SE1.h"
#include "algo.h"

componentRES_SE1::componentRES_SE1(dataRES *d)
:componentLINEAR(){
  if (!d)
    ACE_ERROR("componentRES_SE1::componentRES_SE1 dataRES null");
  mData = (*d);
  mNodePos=mData.nodePos;
  mNodeNeg=mData.nodeNeg;
  mName = mData.name;
  if (ACE_IS_NULL(mData.value))
    ACE_ERROR("Resistor null");
  
  mType = ACE_TYPE_RES;
}
/*stamp:
  
     U=Ri
          U=Vn-Vp
     <-----------
----|------->----|--->
    n       i    p

        | Vp   |  Vn   |
___________________________
KCL(p)  |-1/R  | 1/R   |
___________________________
KCL(n)  |1/R  | -1/R   |
*/
void componentRES_SE1::stamp(){
  //KCLs
  //  int vp,vn;
  //equationKCL *pk;
//   vp=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodePos);
//   vn=algo::spls->getIndexUnknown(ACE_TYPE_V,mData.nodeNeg);
  //pk=;
  ACE_DOUBLE * X = algo::spls->KCL(mData.nodePos)->mCoefs + algo::spls->mDimx;
  ACE_DOUBLE * Zs = algo::spls->KCL(mData.nodePos)->mCoefs + 2*algo::spls->mDimx;
  algo::spls->mNodes[mData.nodePos]->stampV(-1/mData.value,X,Zs);
  algo::spls->mNodes[mData.nodeNeg]->stampV(1/mData.value,X,Zs);

  X = algo::spls->KCL(mData.nodeNeg)->mCoefs + algo::spls->mDimx;
  Zs = algo::spls->KCL(mData.nodeNeg)->mCoefs + 2*algo::spls->mDimx;
  algo::spls->mNodes[mData.nodePos]->stampV(1/mData.value,X,Zs);
  algo::spls->mNodes[mData.nodeNeg]->stampV(-1/mData.value,X,Zs);


  
//   algo::spls->KCL(mData.nodePos)->mCoefs[vp]+=-1/mData.value;
//   algo::spls->KCL(mData.nodePos)->mCoefs[vn]+=1/mData.value;
  
//   algo::spls->KCL(mData.nodeNeg)->mCoefs[vp]+=1/mData.value;
//   algo::spls->KCL(mData.nodeNeg)->mCoefs[vn]+=-1/mData.value;

}
void componentRES_SE1::print(){
  componentLINEAR::print();
  printf("\t value: %f\n",mData.value);
}

componentRES_SE1::~componentRES_SE1(){
  
}
