#include "perform.h"
#include <ngspice.h>
#include "cktdefs.h"
#include "fteext.h"

#include "./spicelib/devices/isrc/isrcdefs.h"
#include "./spicelib/devices/res/resdefs.h"
#include "./spicelib/devices/dio/diodefs.h"
#include "./spicelib/devices/dio/dioext.h"
#include "./spicelib/devices/cap/capext.h"
#include "./spicelib/devices/cap/capdefs.h"
#include "./spicelib/devices/res/resext.h"
#include "./spicelib/devices/ind/inddefs.h"
#include "./spicelib/devices/ind/indext.h"
#include "./spicelib/devices/vsrc/vsrcdefs.h"
#include "./spicelib/devices/vcvs/vcvsdefs.h"
#include "./spicelib/devices/vccs/vccsdefs.h"
#include "./spicelib/devices/asrc/asrcdefs.h"
#include "./spicelib/devices/mos1/mos1defs.h"
//#include "./spicelib/devices/dev.h"
#include "devdefs.h"

#include "./spicelib/devices/dev.h"

static int sType =0;
static char sTypeName[16];
static GENinstance *psInstance =0;

int computeSourcesValues(double time){
  SPICEdev **myDEVices = devices();
  int dev[2];
  int i =0;
  int error=0;
  CKTcircuit *circuit =0;
  dev[0]=46;
  dev[1]=27;
  circuit =(CKTcircuit *) ft_curckt->ci_ckt;

  if (!circuit)
    return 0;
  circuit->CKTtime = time;

  for(i=0;i<2;i++){
    if ( ((*myDEVices[dev[i]]).DEVload != NULL) && (circuit->CKThead[dev[i]] != NULL) ){
      error = (*((*myDEVices[dev[i]]).DEVload))(circuit->CKThead[dev[i]],circuit);
    }
  }
  return 1;
}


int getSourceValue(char *type,void* id,double* value){
  VSRCinstance *hereVSRC;
  ISRCinstance *hereISRC;
  int _type=0;
  _type=INPtypelook(type);
  if (_type == 46){
    hereVSRC = (VSRCinstance *) id;
    *value = hereVSRC->currentValue;
  }
  else if (_type == 27){
    hereISRC = (ISRCinstance *) id;
    *value = hereISRC->currentValue;
  }else{
    return 0;
  }
  return 1;
}

void     fillResistorInfos(void *data, GENinstance *pInstance){
  dataRES *p;
  RESinstance *here;
  if (!data || !psInstance)
    return;
  
  p = (dataRES *)data;
  here = (RESinstance *)psInstance;
  p->name = here->RESname;
  p->nodePos = here->RESposNode;
  p->nodeNeg = here->RESnegNode;
  p->value = here->RESresist;
}
 void     fillInductorInfos(void *data, GENinstance *pInstance){
  dataIND *p;
  INDinstance *here;
  if (!data || !psInstance)
    return;

  p = (dataIND *)data;
  here = (INDinstance *)psInstance;
  p->name = here->INDname;
  p->nodePos = here->INDposNode;
  p->nodeNeg = here->INDnegNode;
  p->value = here->INDinduct;
}
 void     fillVsourceInfos(void *data, GENinstance *pInstance){
  dataVSRC *p;
  VSRCinstance *here;
  if (!data || !psInstance)
    return;
  
  p = (dataVSRC *)data;
  here = (VSRCinstance *)psInstance;
  p->name = here->VSRCname;
  p->nodePos = here->VSRCposNode;
  p->nodeNeg = here->VSRCnegNode;
  p->value = here->VSRCdcValue;
  p->id = (void*)here;
}
void     fillVCVsourceInfos(void *data, GENinstance *pInstance){
  dataVCVS *p;
  VCVSinstance *here;
  if (!data || !psInstance)
    return;
  
  p = (dataVCVS *)data;
  here = (VCVSinstance *)psInstance;
  p->name = here->VCVSname;
  p->nodePos = here->VCVSposNode;
  p->nodeNeg = here->VCVSnegNode;
  p->nodeDriverPos = here->VCVScontPosNode;
  p->nodeDriverNeg = here->VCVScontNegNode;
  p->value = 0;
  p->coef = here->VCVScoeff;
}
void     fillVCCsourceInfos(void *data, GENinstance *pInstance){
  dataVCCS *p;
  VCCSinstance *here;
  if (!data || !psInstance)
    return;
  
  p = (dataVCCS *)data;
  here = (VCCSinstance *)psInstance;
  p->name = here->VCCSname;
  p->nodePos = here->VCCSposNode;
  p->nodeNeg = here->VCCSnegNode;
  p->nodeDriverPos = here->VCCScontPosNode;
  p->nodeDriverNeg = here->VCCScontNegNode;
  p->value = 0;
  p->coef = here->VCCScoeff;
}

 
 void     fillIsourceInfos(void *data, GENinstance *pInstance){
   ISRCinstance * here;
   dataISRC *p;
  if (!data || !psInstance)
    return;
  p = (dataISRC *)data;
  here = (ISRCinstance *)psInstance;
  p->name = here->ISRCname;
  p->nodePos = here->ISRCposNode;
  p->nodeNeg = here->ISRCnegNode;
  p->value = here->ISRCdcValue;
  p->id = (void*)here;
   ;}

 void     fillARBsourceInfos(void *data, GENinstance *pInstance){
   ASRCinstance *here;
   dataARB *p;
  if (!data || !psInstance)
    return;

  p = (dataARB *)data;
  here = (ASRCinstance *)psInstance;
  p->name = here->ASRCname;
  p->nodePos = here->ASRCposNode;
  p->nodeNeg = here->ASRCnegNode;
  p->type = here->ASRCtype;
 }
 void     fillCapacitorInfos(void *data, GENinstance *pInstance){
   CAPinstance *here;
   dataCAP *p;
  if (!data || !psInstance)
    return;

  p = (dataCAP *)data;
  here = (CAPinstance *)psInstance;
  p->name = here->CAPname;
  p->nodePos = here->CAPposNode;
  p->nodeNeg = here->CAPnegNode;
  p->value = here->CAPcapac;
 }
 void     fillDiodeInfos(void *data, GENinstance *pInstance){
   DIOinstance *here;
   dataDIO *p;
  if (!data || !psInstance)
    return;
  p = (dataDIO *)data;
  here = (DIOinstance *)psInstance;
  p->name = here->DIOname;
  p->nodePos = here->DIOposNode;
  p->nodeNeg = here->DIOnegNode;
  p->area = here->DIOarea;
 }
void bidon(MOS1instance *here){
  printf("%d\n",here->MOS1mode);
}
 void fillMos1Infos(void *data, GENinstance *pInstance){
   MOS1instance *here;
   dataMOS1 *p;
  if (!data || !psInstance)
    return;
  p = (dataMOS1 *)data;
  here = (MOS1instance *)psInstance;
  bidon(here);
  p->mode = here->sMOS1modPtr->MOS1type;
  p->name = here->MOS1name;
  p->drain = here->MOS1dNode;
  p->gate = here->MOS1gNode;
  p->source = here->MOS1sNode;
  p->w =  here->MOS1w;
  p->k =  here->sMOS1modPtr->MOS1transconductance;
  p->vt = here->sMOS1modPtr->MOS1vt0;

  
 }

void MEperform(){
    CKTcircuit *circuit =0;
    circuit =(CKTcircuit *) ft_curckt->ci_ckt;
    
}

void MEprint(){
  CKTcircuit *circuit =0;
  int type =0;
  GENmodel *pModel =0;
  GENinstance *pInstance =0;
  ISRCinstance *pISRC =0;
  VSRCinstance * pVSRC =0;
  double val =0;
  CKTnode * node;

  int nodeNumber=0;
  if (!ft_curckt || !ft_curckt->ci_ckt){
    printf("Pas de circuit\n");
    return;
  }
  
  printf("******DEBUT DE LA LECTURE DE LA TABLE******\n");
  circuit =(CKTcircuit *) ft_curckt->ci_ckt;
  
  printf("******LES %d NODES :******\n",getNbElementsOfType("Node"));

  node = circuit->CKTnodes;
  while (node){
    nodeNumber = node->number;
    printf("node number: %d\n \t name: %s\n \t type: %d \n",nodeNumber,node->name,node->type);
    node = node->next;
  }
  printf("-Les sources d'intensite:\n");
  printf("-------------------------\n");
  type = INPtypelook("Isource");
  pModel =  circuit->CKThead[type];
  if (pModel){
    pInstance = pModel->GENinstances;
    while (pInstance != 0){
      pISRC = (ISRCinstance *) pInstance;
      val = pISRC->ISRCdcValue;
      printf("une source entre les noeuds neg %d et pos %d, de valeur %f A.\n",pISRC->ISRCnegNode,pISRC->ISRCposNode,val);
      pInstance = pInstance->GENnextInstance;
    }
  }
  
  printf("-Les capacitees:\n");
  printf("---------------\n");
  type = INPtypelook("Capacitor");
  pModel =  circuit->CKThead[type];
  CAPsPrint(pModel,circuit);

  printf("les diodes: \n");
  printf("-----------");
  type = INPtypelook("Diode");
  pModel =  circuit->CKThead[type];
  DIOsPrint(pModel, circuit);

  printf("les resistances:\n");
  printf("----------------");
  type = INPtypelook("Resistor");
  pModel =  circuit->CKThead[type];
  RESsPrint(pModel, circuit);

  printf("independent voltage sources\n");
  printf("----------------------------\n");
  type = INPtypelook("Vsource");
  pModel = circuit->CKThead[type];
  if (pModel){
    pInstance = pModel->GENinstances;
    while (pInstance != 0){
      pVSRC = (VSRCinstance *) pInstance;
      val = pVSRC->VSRCdcValue;
      printf("name : %s\n entre les noeuds neg %d et pos %d.\n valeur : %f Volt \n",pVSRC->VSRCname,  pVSRC->VSRCnegNode,pVSRC->VSRCposNode,val);
      pInstance = pInstance->GENnextInstance;
    }
  }

  printf("Les inductors:\n");
  printf("--------------\n");
  type = INPtypelook("Inductor");
  INDsPrint(circuit->CKThead[type],circuit);
   
  printf("******FIN DE LA LECTURE DE LA TABLE******\n");
}
int initComponentList(char *type){
  GENmodel *pModel =0;
  CKTcircuit *circuit =0;
  circuit =(CKTcircuit *) ft_curckt->ci_ckt;

  if (!circuit)
    return -1;
  strcpy(sTypeName,type);
  sType = INPtypelook(type);
  pModel =  circuit->CKThead[sType];
  if (!pModel)
    return 0;
  psInstance = pModel->GENinstances;
  if (!psInstance)
    return 0;
  return 1;
}
int nextComponent(void * data){
  if (!psInstance)
    return 0;
  switch(sType){
  case 38:
    fillResistorInfos(data,psInstance);
    break;
  case 25:
    fillInductorInfos(data,psInstance);
    break;
  case 46:
    fillVsourceInfos(data,psInstance);
    break;
  case 27:
    fillIsourceInfos(data,psInstance);
    break;
  case 21:
    fillDiodeInfos(data,psInstance);
    break;
  case 33:
    fillMos1Infos(data,psInstance);
    break;
  case 16:
    fillCapacitorInfos(data,psInstance);
    break;
  case 45:
    fillVCVsourceInfos(data,psInstance);
    break;
  case 1:
    fillARBsourceInfos(data,psInstance);
    break;
  case 44:
    fillVCCsourceInfos(data,psInstance);
    break;
  default :
    printf("ERROR parser/src/perform.c : unknown type\n");

    
    return 0;
  }
  psInstance = psInstance->GENnextInstance;
  return 1;
}
int getNbElementsOfType(char* type){
  
  CKTcircuit *circuit =0;
  circuit =(CKTcircuit *) ft_curckt->ci_ckt;
  if (!strcmp(type, "Node")){
    if (circuit && circuit->CKTlastNode)
      return circuit->CKTlastNode->number +1;
    else{
      printf("parser : getNbElementsOfType no node\n");
      return -1;
    }
  }
  else if(!strcmp(type, "Capacitor")){
    return sNbCap;
  }
  printf("parser : getNbElementsOfType not implemented %s \n",type);
  return -1;
}

