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

/*! \file acetime.h

*/
/************************************************************************
  			acetime.h 

**************************************************************************/
#include <stdio.h>
#include "sys/time.h"
#include <iostream>
using namespace std;


#ifndef ACETIME_H
#define ACETIME_H
// Class aceTime
// 
// 
class aceTime {
public:
  aceTime();
  void start();
  void stop();
  void setName(char *Name);
  void print(ostream& os = cout);
  virtual ~aceTime();
  
protected:
private:
  long mCall;
  timeval mStart;
  long mCumul;
  char mName[128];
  bool mIsRunning;
    

};
#endif //ACETIME_H

