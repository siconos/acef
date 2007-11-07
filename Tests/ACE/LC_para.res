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
	c1_U1_0	l1_I0_1
Zs
	V0	V1
Zns

---------------------------------------
equation	c1_U1_0'	l1_I0_1'	c1_U1_0	l1_I0_1	V0	V1
KCL1*0		-0.01	0	0	1	0	0	1
IND*1		0	0.02	0	0	1	-1	0
KCL0		0	0	0	-1	0	0	-1
TEN		0	0	-1	0	1	-1	0
Ax'=Bx+CZs+DZns+s
-----------------
A:
[2,2]
	-0.01	0
	0	0.02
B:
[2,2]
	0	1
	0	0
C:
[2,1]
	0
	-1
D:
s:
[2,1]
	1
	0
inv A:
-----
[2,2]
	-100	0
	0	50
system x'=A1x+A1zs+A1zns+s
A1x:
[2,2]
	0	-100
	0	0
A1zs:
[2,1]
	0
	-50
A1zns:
A1s:
[2,1]
	-100
	0
print cap i coefs c1:
	0.000000	0.000000	0.000000	-1.000000	0.000000	0.000000	-1.000000
ACE MESSAGE :final equation ;
--->linearSystem with 3 equations whose 2 dynamic equations.
x
	c1_U1_0	l1_I0_1
Zs
	V0	V1
Zns

---------------------------------------
equation	c1_U1_0'	l1_I0_1'	c1_U1_0	l1_I0_1	V0	V1
KCL1*0		-0.01	0	0	1	0	0	1
IND*1		0	0.02	0	0	1	-1	0
KCL0		0	0	0	0	0	0	0
TEN		0	0	-1	0	1	-1	0
0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
B1x:
[1,2]
	-1	0
B1zs:
[1,1]
	-1
B1zns:
B1s:
[1,1]
	0
Zns = C1x*x + C1s*Zs + C1l*lamdba + C1s
C1x:
C1zs:
C1l:
C1s:
Y = D1x*x + D1s*Zs + D1ns*Zns + D1l*lambda +D1s
D1x:
D1zs:
D1zns:
D1l:
D1s:
R=A1zns*C1l
x'=A2x*x + A2zs*Zs + R*lambda+A2s
0=B2x*x + B2zs*Zs + B2l*lambda + B2s
Y=D2x*x + D2zs*Zs + D2l*lambda + D2s
R:
A2x:
[2,2]
	0	-100
	0	0
A2zs:
[2,1]
	0
	-50
A2s:
[2,1]
	-100
	0
B2x:
[1,2]
	-1	0
B2zs:
[1,1]
	-1
B2l:
[1,0]

B2s:
[1,1]
	0
D2x:
D2zs:
D2l:
D2s:
