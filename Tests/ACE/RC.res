rep1
rep2
readfile
main rep6 fopen : RC.cir

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
 entre les noeuds neg 0 et pos 1.
 valeur : 3.000000 Volt 
Les inductors:
--------------
INDUCTORS----------
******FIN DE LA LECTURE DE LA TABLE******
print edge in:
print edge out:
ACE 0 Inductors:
ACE 1 Capacitors:
component type : Capacitor 
	 Name : c1
	nodePos, nodeNeg : 2 0
ACE 1 Resistors:
component type : Resistor 
	 Name : r1
	nodePos, nodeNeg : 1 2
	 value: 100.000000
ACE 0 Isource:
ACE 1 Vsource:
component type : Vsource 
	 Name : vs
	nodePos, nodeNeg : 1 0
ACE 0 Diode:
--->linearSystem with 4 equations whose 1 dynamic equations.
x
	c1_U2_0
Zs
	V0	V1	V2	vs_I0_1
Zns

---------------------------------------
equation	c1_U2_0'	c1_U2_0	V0	V1	V2	vs_I0_1
KCL2*0		-0.02	0	0	0.01	-0.01	0	0
KCL0		0	0	0	0	0	-1	0
KCL1		0	0	0	-0.01	0.01	1	0
VD		0	0	-1	1	0	0	-3
TEN		0	-1	1	0	-1	0	0
Ax'=Bx+CZs+DZns+s
-----------------
A:
[1,1]
	-0.02
B:
[1,1]
	0
C:
[1,4]
	0	0.01	-0.01	0
D:
s:
[1,1]
	0
inv A:
-----
[1,1]
	-50
system x'=A1x+A1zs+A1zns+s
A1x:
[1,1]
	0
A1zs:
[1,4]
	0	-0.5	0.5	0
A1zns:
A1s:
[1,1]
	0
print cap i coefs c1:
	0.000000	0.000000	0.000000	-0.010000	0.010000	0.000000	0.000000
ACE_MESSAGE :final equation ;
--->linearSystem with 4 equations whose 1 dynamic equations.
x
	c1_U2_0
Zs
	V0	V1	V2	vs_I0_1
Zns

---------------------------------------
equation	c1_U2_0'	c1_U2_0	V0	V1	V2	vs_I0_1
KCL2*0		-0.02	0	0	0.01	-0.01	0	0
KCL0		0	0	0	0.01	-0.01	-1	0
KCL1		0	0	0	-0.01	0.01	1	0
VD		0	0	-1	1	0	0	-3
TEN		0	-1	1	0	-1	0	0
0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
B1x:
[3,1]
	0
	0
	-1
B1zs:
[3,4]
	0	-0.01	0.01	1
	-1	1	0	0
	1	0	-1	0
B1zns:
B1s:
[3,1]
	0
	-3
	0
