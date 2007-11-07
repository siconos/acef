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
--->linearSystem with 3 equations whose 1 dynamic equations.
x
	l1_I0_2
Zs
	V0	V1	V2
Zns

---------------------------------------
equation	l1_I0_2'	l1_I0_2	V0	V1	V2
IND*0		0.02	0	1	0	-1	0
KCL0		0	-1	0	0	0	-1
KCL1		0	0	0	-0.1	0.1	1
KCL2		0	1	0	0.1	-0.1	0
Ax'=Bx+CZs+DZns+s
-----------------
A:
[1,1]
	0.02
B:
[1,1]
	0
C:
[1,2]
	0	-1
D:
s:
[1,1]
	0
inv A:
-----
[1,1]
	50
system x'=A1x+A1zs+A1zns+s
A1x:
[1,1]
	0
A1zs:
[1,2]
	0	-50
A1zns:
A1s:
[1,1]
	0
ACE MESSAGE :final equation ;
--->linearSystem with 3 equations whose 1 dynamic equations.
x
	l1_I0_2
Zs
	V0	V1	V2
Zns

---------------------------------------
equation	l1_I0_2'	l1_I0_2	V0	V1	V2
IND*0		0.02	0	1	0	-1	0
KCL0		0	-1	0	0	0	-1
KCL1		0	0	0	-0.1	0.1	1
KCL2		0	1	0	0.1	-0.1	0
0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
B1x:
[2,1]
	0
	1
B1zs:
[2,2]
	-0.1	0.1
	0.1	-0.1
B1zns:
B1s:
[2,1]
	1
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
[1,1]
	0
A2zs:
[1,2]
	0	-50
A2s:
[1,1]
	0
B2x:
[2,1]
	0
	1
B2zs:
[2,2]
	-0.1	0.1
	0.1	-0.1
B2l:
[2,0]


B2s:
[2,1]
	1
	0
D2x:
D2zs:
D2l:
D2s:
