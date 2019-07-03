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

/*! \file unknown.cpp

*/
/************************************************************************
  			unknown.cpp 
**************************************************************************/
#include "unknown.h"


unknown::unknown(int type, component *c)
{
  mNode =-1;
  mIndex=-1;
  mDynIndex=-1;
  mIndexInVector = -1;
  mType = type;
  mComponent = c;
  buildName();
}
unknown::unknown(int type, int node)
{
  mNode =node;
  mIndex=-1;
  mDynIndex=-1;
  mType = type;
  mComponent = 0;
  buildName();
}
/*useful only for display and debug*/
void unknown::buildName(){
  char type[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,type);
  if (mType == ACE_TYPE_V){
    sprintf(mName,"V%d",mNode);
  }else{
    if (mComponent){
      if ((mComponent->mType == ACE_TYPE_COMP || mComponent->mType == ACE_TYPE_RELAY)){
	if (mType == ACE_TYPE_I)
	  sprintf(mName,"I_%s_S",mComponent->mName);
	else
	  sprintf(mName,"%s_%s_S",type,mComponent->mName);
      }else{
	if (mType == ACE_TYPE_I)
	  sprintf(mName,"I_%s_%d_%d",mComponent->mName,mComponent->mNodeNeg,mComponent->mNodePos);
	else
	  sprintf(mName,"%s_%s_%d_%d",type,mComponent->mName,mComponent->mNodeNeg,mComponent->mNodePos);
      }
    }else{
      sprintf(mName,"%s",type);
    }
  }
}
void unknown::print(){
  printf("\t%s",mName);  
}

void unknown::printdev()
{
  print();
  printf("'");
}
