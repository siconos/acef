rep1
rep2
readfile
main rep6 fopen : RLC.cir

Circuit: LINEAR_DC_CKT.CIR - SIMPLE CIRCUIT FOR NODAL ANALYSIS

readfile end
******DEBUT DE LA LECTURE DE LA TABLE******
******LES 4 NODES :******
node number: 0
 	 name: 0
 	 type: 3 
node number: 1
 	 name: 1
 	 type: 3 
node number: 2
 	 name: 2
 	 type: 3 
node number: 3
 	 name: 3
 	 type: 3 
-Les sources d'intensite:
-------------------------
-Les capacitees:
---------------
CAPACITORS-----------------
Model name:C
    Instance name:c1
      Positive, negative nodes: 2, 3
      Capacitance: 2.000000e-02(specified)
    CAPsenParmNo:0
les diodes: 
-----------DIOS-----------------
les resistances:
----------------RESISTORS-----------------
Model name:R
    Instance name:r1
      Positive, negative nodes name: 1, 2
      Positive, negative nodes: 1, 2
  Multiplier: 0 (default)
      Resistance: 100.000000 (specified)
    RESsenParmNo:0
independent voltage sources
----------------------------
name : vs
 entre les noeuds neg 1 et pos 0.
 valeur : 3.000000 Volt 
Les inductors:
--------------
INDUCTORS----------
Model name:L
    Instance name:l1
      Positive, negative nodes: 3, 0
      Branch Equation: 0
      Inductance: 0.03 (specified)
    INDsenParmNo:0
******FIN DE LA LECTURE DE LA TABLE******
print edge in:
2  3
print edge out:
ACE 1 Inductors:
component type : Inductor 
	 Name : l1
	nodePos, nodeNeg : 3 0
ACE 1 Capacitors:
component type : Capacitor 
	 Name : c1
	nodePos, nodeNeg : 2 3
ACE 1 Resistors:
component type : Resistor 
	 Name : r1
	nodePos, nodeNeg : 1 2
	 value: 100.000000
ACE 0 Isource:
ACE 1 Vsource:
component type : Vsource 
	 Name : vs
	nodePos, nodeNeg : 0 1
ACE 0 Diode:
--->linearSystem with 6 equations whose 2 dynamic equations.
x
	c1_U2_3	l1_I0_3
Zs
	V0	V1	V2	V3	vs_I1_0
Zns

---------------------------------------
equation	c1_U2_3'	l1_I0_3'	c1_U2_3	l1_I0_3	V0	V1	V2	V3	vs_I1_0
KCL2*0		-0.02	0	0	0	0	0.01	-0.01	0	0	0
IND*1		0	0.03	0	0	1	0	0	-1	0	0
KCL0		0	0	0	-1	0	0	0	0	1	0
KCL1		0	0	0	0	0	-0.01	0.01	0	-1	0
KCL3		0	0	0	1	0	0	0	0	0	0
VD		0	0	0	0	1	-1	0	0	0	-3
TEN		0	0	-1	0	0	0	-1	1	0	0
Ax'=Bx+CZs+DZns+s
-----------------
A:
[2,2]
	-0.02	0
	0	0.03
B:
[2,2]
	0	0
	0	0
C:
[2,5]
	0	0.01	-0.01	0	0
	1	0	0	-1	0
D:
s:
[2,1]
	0
	0
inv A:
-----
[2,2]
	-50	0
	0	33.3333
system x'=A1x+A1zs+A1zns+s
A1x:
[2,2]
	0	0
	0	0
A1zs:
[2,5]
	0	-0.5	0.5	0	0
	33.3333	0	0	-33.3333	0
A1zns:
A1s:
[2,1]
	0
	0
print cap i coefs c1:
	0.000000	0.000000	0.000000	0.000000	0.000000	-0.010000	0.010000	0.000000	0.000000	0.000000
ACE_MESSAGE :final equation ;
--->linearSystem with 6 equations whose 2 dynamic equations.
x
	c1_U2_3	l1_I0_3
Zs
	V0	V1	V2	V3	vs_I1_0
Zns

---------------------------------------
equation	c1_U2_3'	l1_I0_3'	c1_U2_3	l1_I0_3	V0	V1	V2	V3	vs_I1_0
KCL2*0		-0.02	0	0	0	0	0.01	-0.01	0	0	0
IND*1		0	0.03	0	0	1	0	0	-1	0	0
KCL0		0	0	0	-1	0	0	0	0	1	0
KCL1		0	0	0	0	0	-0.01	0.01	0	-1	0
KCL3		0	0	0	1	0	0.01	-0.01	0	0	0
VD		0	0	0	0	1	-1	0	0	0	-3
TEN		0	0	-1	0	0	0	-1	1	0	0
0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
B1x:
[4,2]
	0	0
	0	1
	0	0
	-1	0
B1zs:
[4,5]
	0	-0.01	0.01	0	-1
	0	0.01	-0.01	0	0
	1	-1	0	0	0
	0	0	-1	1	0
B1zns:
B1s:
[4,1]
	0
	0
	-3
	0
