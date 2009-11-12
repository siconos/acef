/************************************************************************
  			componentres_SE1.h 

**************************************************************************/

#ifndef COMPONENTRES_SE1_H
#define COMPONENTRES_SE1_H
#include "ace.h"

#include "componentlinear.h"

/** The componentres class.
*/


/**
 *
 *It corresponds to a resistor component.\n
 *
**************************************************************************/


// Class componentRES_SE1
// 
// 
class componentRES_SE1 : public componentLINEAR {

public:
  dataRES mData;
  componentRES_SE1(dataRES *d);
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
  virtual ~componentRES_SE1();

protected:
private:
};
#endif //COMPONENTRES_SE1_H

