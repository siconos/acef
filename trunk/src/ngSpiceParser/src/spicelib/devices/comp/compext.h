/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 2008 O. Bonnefon
**********/

extern int COMPask(CKTcircuit*,GENinstance*,int,IFvalue*,IFvalue*);
extern int COMPdelete(GENmodel*,IFuid,GENinstance**);
extern void COMPdestroy(GENmodel**);
extern int COMPload(GENmodel*,CKTcircuit*);
extern int COMPacload(GENmodel*,CKTcircuit*);
extern int COMPmodAsk(CKTcircuit*,GENmodel*,int,IFvalue*);
extern int COMPmDelete(GENmodel**,IFuid,GENmodel*);
extern int COMPmParam(int,IFvalue*,GENmodel*);
extern int COMPparam(int,IFvalue*,GENinstance*,IFvalue*);
extern int COMPpzLoad(GENmodel*,CKTcircuit*,SPcomplex*);
extern int COMPsAcLoad(GENmodel*,CKTcircuit*);
extern int COMPsLoad(GENmodel*,CKTcircuit*);
extern int COMPsSetup(SENstruct*,GENmodel*);
extern void COMPsPrint(GENmodel*,CKTcircuit*);
extern int COMPsetup(SMPmatrix*,GENmodel*,CKTcircuit*,int*);
extern int COMPtemp(GENmodel*,CKTcircuit*);
extern int COMPnoise(int,int,GENmodel*,CKTcircuit*,Ndata*,double*);
