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

static int sType =0;
static char sTypeName[16];
static GENinstance *psInstance =0;

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
   ;}
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
   ;}
 

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
      printf("une source entre les noeuds %d et %d, de valeur %f A.\n",pInstance->GENnode1,pInstance->GENnode2,val);
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
      printf("name : %s\n entre les noeuds %d et %d.\n valeur : %f Volt \n",pVSRC->VSRCname,  pInstance->GENnode1,pInstance->GENnode2,val);
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
    case 18:
      fillVsourceInfos(data,psInstance);
      break;
    case 27:
      fillIsourceInfos(data,psInstance);
      break;
    case 21:
      fillDiodeInfos(data,psInstance);
      break;
    case 16:
      fillCapacitorInfos(data,psInstance);
      break;
    default :
      /*error*/;
    return -2;
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

