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

/*! \file acetime.cpp

*/
/************************************************************************
  			acetime.cpp 
**************************************************************************/


#include "acetime.h"
#include "string.h"

#define ACE_WITH_TIMER

aceTime::~aceTime(){
}
aceTime::aceTime(){
#ifdef ACE_WITH_TIMER
  mIsRunning=false;
  mCumul=0;
  mCall=0;
  strcpy(mName,"NoName");
#endif
}

void aceTime::start(){
#ifdef ACE_WITH_TIMER
  gettimeofday(&mStart,NULL);
  mCall++;
  mIsRunning=true;
#endif
}

void aceTime::stop(){
#ifdef ACE_WITH_TIMER
  if (mIsRunning){
    timeval aux;
    gettimeofday(&aux,NULL);
    mCumul += (aux.tv_sec - mStart.tv_sec)*1000000 +(aux.tv_usec - mStart.tv_usec) ;
    mIsRunning = false;
  }
#endif
}
void aceTime::print(ostream& os){
#ifdef ACE_WITH_TIMER
  os << mName<<" : "<<mCumul<<" usec. "<<mCall<<" calls.";
  if (mCall)
    os << "ie "<<mCumul/mCall<<" per call.";
  os <<endl;
#endif
}
void aceTime::setName(char *Name){
#ifdef ACE_WITH_TIMER
  strcpy(mName,Name);
#endif
}


