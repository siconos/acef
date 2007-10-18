/************************************************************************
  			unknown.cpp 
**************************************************************************/
#include "unknown.h"


unknown::unknown(int type, component *c)
{
  mNode =-1;
  mIndex=-1;
  mDynIndex=-1;
  mType = type;
  mComponent = c;
}
unknown::unknown(int type, int node)
{
  mNode =node;
  mIndex=-1;
  mDynIndex=-1;
  mType = type;
  mComponent = 0;
}
void unknown::print(){
  char type[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,type);
  if (mType == ACE_TYPE_V){
    printf("\t%s%d",type,mNode);
  }else{
    if (mComponent){
      if (mType == ACE_TYPE_I)
	printf("\t%s%d_%d",type,mComponent->mNodeNeg,mComponent->mNodePos);
      else
	printf("\t%s%d_%d",type,mComponent->mNodePos,mComponent->mNodeNeg);
    }else{
      printf("\t%s",type);
    }
  }
  
}

