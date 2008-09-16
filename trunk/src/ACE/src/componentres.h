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

