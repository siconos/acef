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

/*! \file componentres.h

*/
/************************************************************************
  			componentres.h 

**************************************************************************/

#ifndef COMPONENTRES_H
#define COMPONENTRES_H
#include "ace.h"

#include "componentlinear.h"

/** The componentres class.
*/


/**
 *
 *It corresponds to a resistor component.\n
 *
**************************************************************************/


// Class componentRES
// 
// 
class componentRES : public componentLINEAR {

public:
  dataRES mData;
  componentRES(dataRES *d);
  /**
     This method fills the tableau equations\n
 <TABLE >
          <CAPTION>Stamp method</CAPTION>
          <TR>
             <TH> </TH>
             <TH>Vp</TH>
             <TH>Vn</TH>
          </TR>
          <TR>
             <TD>KCL(p)</TD>
             <TD>-1/R</TD>
             <TD>1/R</TD>
          </TR>
          <TR>
             <TD>KCL(n)</TD>
             <TD>1/R</TD>
             <TD>-1/R</TD>
          </TR>
        </TABLE>
   */

  virtual void stamp();
  virtual void print();
  virtual ~componentRES();

protected:
private:
};
#endif //COMPONENTRES_H

