rep1
rep2
readfile
main rep6 fopen : LR.cir

Circuit: LINEAR_DC_CKT.CIR - SIMPLE CIRCUIT FOR NODAL ANALYSIS

readfile end
******DEBUT DE LA LECTURE DE LA TABLE******
******LES 3 NODES :******
node number: 0
 	 name: 0
 	 type: 3 
node number: 1
 	 name: 1
 	 type: 3 
node number: 2
 	 name: 2
 	 type: 3 
-Les sources d'intensite:
-------------------------
une source entre les noeuds neg 0 et pos 1, de valeur 1.000000 A.
-Les capacitees:
---------------
CAPACITORS-----------------
les diodes: 
-----------DIOS-----------------
les resistances:
----------------RESISTORS-----------------
Model name:R
    Instance name:r1
      Positive, negative nodes name: 1, 2
      Positive, negative nodes: 1, 2
  Multiplier: 0 (default)
      Resistance: 10.000000 (specified)
    RESsenParmNo:0
independent voltage sources
----------------------------
Les inductors:
--------------
INDUCTORS----------
Model name:L
    Instance name:l1
      Positive, negative nodes: 2, 0
      Branch Equation: 0
      Inductance: 0.02 (specified)
    INDsenParmNo:0
******FIN DE LA LECTURE DE LA TABLE******
x
	________I0_2
Zs
	________V0	________V1	________V2
Zns

---------------------------------------
equation	________I0_2	________I0_2	________V0	________V1	________V2
KCL0		0.000000000000	-1.000000000000	0.000000000000	0.000000000000	0.000000000000	-1.000000000000
KCL1		0.000000000000	0.000000000000	0.000000000000	-0.100000000000	0.100000000000	1.000000000000
KCL2		0.000000000000	1.000000000000	0.000000000000	0.100000000000	-0.100000000000	0.000000000000
IND*		0.020000000000	0.000000000000	1.000000000000	0.000000000000	-1.000000000000	0.000000000000
ACE 1 Inductors:
component type : Inductor 
	 Name : l1
	nodePos, nodeNeg : 2 0
ACE 0 Capacitors:
ACE 1 Resistors:
component type : Resistor 
	 Name : r1
	nodePos, nodeNeg : 1 2
	 value: 10.000000
ACE 1 Isource:
component type : Isource 
	 Name : is
	nodePos, nodeNeg : 1 0
	 value: 1.000000
ACE 0 Vsource:
ACE 0 Diode:
