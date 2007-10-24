rep1
rep2
readfile
main rep6 fopen : LC_para.cir

Circuit: LINEAR_DC_CKT.CIR - SIMPLE CIRCUIT FOR NODAL ANALYSIS

readfile end
******DEBUT DE LA LECTURE DE LA TABLE******
******LES 2 NODES :******
node number: 0
 	 name: 0
 	 type: 3 
node number: 1
 	 name: 1
 	 type: 3 
-Les sources d'intensite:
-------------------------
une source entre les noeuds neg 0 et pos 1, de valeur 1.000000 A.
-Les capacitees:
---------------
CAPACITORS-----------------
Model name:C
    Instance name:c1
      Positive, negative nodes: 1, 0
      Capacitance: 1.000000e-02(specified)
    CAPsenParmNo:0
les diodes: 
-----------DIOS-----------------
les resistances:
----------------RESISTORS-----------------
independent voltage sources
----------------------------
Les inductors:
--------------
INDUCTORS----------
Model name:L
    Instance name:l1
      Positive, negative nodes: 1, 0
      Branch Equation: 0
      Inductance: 0.02 (specified)
    INDsenParmNo:0
******FIN DE LA LECTURE DE LA TABLE******
print edge in:
print edge out:
ACE 1 Inductors:
component type : Inductor 
	 Name : l1
	nodePos, nodeNeg : 1 0
ACE 1 Capacitors:
component type : Capacitor 
	 Name : c1
	nodePos, nodeNeg : 1 0
ACE 0 Resistors:
ACE 1 Isource:
component type : Isource 
	 Name : is
	nodePos, nodeNeg : 1 0
	 value: 1.000000
ACE 0 Vsource:
ACE 0 Diode:
--->linearSystem with 3 equations whose 2 dynamic equations.
x
	________U1_0	________I0_1
Zs
	________V0	________V1
Zns

---------------------------------------
equation	________U1_0	________I0_1	________U1_0	________I0_1	________V0	________V1
KCL1*0		-0.010000000000	0.000000000000	0.000000000000	1.000000000000	0.000000000000	0.000000000000	1.000000000000
IND*1		0.000000000000	0.020000000000	0.000000000000	0.000000000000	1.000000000000	-1.000000000000	0.000000000000
KCL0		0.000000000000	0.000000000000	0.000000000000	-1.000000000000	0.000000000000	0.000000000000	-1.000000000000
TEN		0.000000000000	0.000000000000	1.000000000000	0.000000000000	1.000000000000	-1.000000000000	0.000000000000
A:
[2,2]
	-0.010000	0.000000
	0.000000	0.020000
inv A:
[2,2]
	-100.000000	0.000000
	0.000000	50.000000
B:
[2,2]
	0.000000	1.000000
	0.000000	0.000000
C:
[2,2]
	0.000000	0.000000
	1.000000	-1.000000
D:
s:
[2,1]
	1.000000
	0.000000
system x'=A1x+A1zs+A1zns+s
A1x:
[2,2]
	0.000000	-100.000000
	0.000000	0.000000
A1zs:
[2,2]
	0.000000	0.000000
	50.000000	-50.000000
A1zns:
A1s:
[2,1]
	-100.000000
	0.000000
print cap i coefs c1:
	0.000000	0.000000	0.000000	-100.000000	0.000000	0.000000	-100.000000
