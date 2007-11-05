rep1
rep2
readfile
main rep6 fopen : LC.cir

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
Model name:C
    Instance name:c1
      Positive, negative nodes: 2, 0
      Capacitance: 2.000000e-02(specified)
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
      Positive, negative nodes: 1, 2
      Branch Equation: 0
      Inductance: 0.01 (specified)
    INDsenParmNo:0
******FIN DE LA LECTURE DE LA TABLE******
print edge in:
print edge out:
ACE 1 Inductors:
component type : Inductor 
	 Name : l1
	nodePos, nodeNeg : 1 2
ACE 1 Capacitors:
component type : Capacitor 
	 Name : c1
	nodePos, nodeNeg : 2 0
ACE 0 Resistors:
ACE 1 Isource:
component type : Isource 
	 Name : is
	nodePos, nodeNeg : 1 0
	 value: 1.000000
ACE 0 Vsource:
ACE 0 Diode:
--->linearSystem with 4 equations whose 2 dynamic equations.
x
	c1_U2_0	l1_I2_1
Zs
	V0	V1	V2
Zns

---------------------------------------
equation	c1_U2_0'	l1_I2_1'	c1_U2_0	l1_I2_1	V0	V1	V2
KCL2*0		-0.02	0	0	-1	0	0	0	0
IND*1		0	0.01	0	0	0	-1	1	0
KCL0		0	0	0	0	0	0	0	-1
KCL1		0	0	0	1	0	0	0	1
TEN		0	0	-1	0	1	0	-1	0
Ax'=Bx+CZs+DZns+s
-----------------
A:
[2,2]
	-0.02	0
	0	0.01
B:
[2,2]
	0	-1
	0	0
C:
[2,3]
	0	0	0
	0	-1	1
D:
s:
[2,1]
	0
	0
inv A:
-----
[2,2]
	-50	0
	0	100
system x'=A1x+A1zs+A1zns+s
A1x:
[2,2]
	0	50
	0	0
A1zs:
[2,3]
	0	0	0
	0	-100	100
A1zns:
A1s:
[2,1]
	0
	0
print cap i coefs c1:
	0.000000	0.000000	0.000000	1.000000	0.000000	0.000000	0.000000	0.000000
ACE_MESSAGE :final equation ;
--->linearSystem with 4 equations whose 2 dynamic equations.
x
	c1_U2_0	l1_I2_1
Zs
	V0	V1	V2
Zns

---------------------------------------
equation	c1_U2_0'	l1_I2_1'	c1_U2_0	l1_I2_1	V0	V1	V2
KCL2*0		-0.02	0	0	-1	0	0	0	0
IND*1		0	0.01	0	0	0	-1	1	0
KCL0		0	0	0	-1	0	0	0	-1
KCL1		0	0	0	1	0	0	0	1
TEN		0	0	-1	0	1	0	-1	0
0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
B1x:
[2,2]
	0	1
	-1	0
B1zs:
[2,3]
	0	0	0
	1	0	-1
B1zns:
B1s:
[2,1]
	1
	0
