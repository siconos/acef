rep1
rep2
readfile
main rep6 fopen : RLC_para.cir

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
      Positive, negative nodes: 0, 1
      Capacitance: 1.000000e-03(specified)
    CAPsenParmNo:0
les diodes: 
-----------DIOS-----------------
les resistances:
----------------RESISTORS-----------------
Model name:R
    Instance name:r1
      Positive, negative nodes name: 0, 1
      Positive, negative nodes: 0, 1
  Multiplier: 0 (default)
      Resistance: 1000.000000 (specified)
    RESsenParmNo:0
independent voltage sources
----------------------------
name : vs
 entre les noeuds neg 1 et pos 0.
 valeur : 2.000000 Volt 
Les inductors:
--------------
INDUCTORS----------
Model name:L
    Instance name:l1
      Positive, negative nodes: 1, 0
      Branch Equation: 0
      Inductance: 0.002 (specified)
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
	nodePos, nodeNeg : 0 1
ACE 1 Resistors:
component type : Resistor 
	 Name : r1
	nodePos, nodeNeg : 0 1
	 value: 1000.000000
ACE 1 Isource:
component type : Isource 
	 Name : is
	nodePos, nodeNeg : 1 0
	 value: 1.000000
ACE 1 Vsource:
component type : Vsource 
	 Name : vs
	nodePos, nodeNeg : 0 1
ACE 0 Diode:
--->linearSystem with 4 equations whose 2 dynamic equations.
x
	________U0_1	________I0_1
Zs
	________V0	________V1	________I1_0
Zns

---------------------------------------
equation	________U0_1	________I0_1	________U0_1	________I0_1	________V0	________V1	________I1_0
KCL1*0		0.001000000000	0.000000000000	0.000000000000	1.000000000000	0.001000000000	-0.001000000000	-1.000000000000	1.000000000000
IND*1		0.000000000000	0.002000000000	0.000000000000	0.000000000000	1.000000000000	-1.000000000000	0.000000000000	0.000000000000
KCL0		0.000000000000	0.000000000000	0.000000000000	-1.000000000000	-0.001000000000	0.001000000000	1.000000000000	-1.000000000000
VD		0.000000000000	0.000000000000	0.000000000000	0.000000000000	1.000000000000	-1.000000000000	0.000000000000	-2.000000000000
TEN		0.000000000000	0.000000000000	1.000000000000	0.000000000000	-1.000000000000	1.000000000000	0.000000000000	0.000000000000
A:
[2,2]
	0.001000	0.000000
	0.000000	0.002000
inv A:
[2,2]
	1000.000000	0.000000
	0.000000	500.000000
B:
[2,2]
	0.000000	1.000000
	0.000000	0.000000
C:
[2,3]
	0.001000	-0.001000	-1.000000
	1.000000	-1.000000	0.000000
D:
s:
[2,1]
	1.000000
	0.000000
system x'=A1x+A1zs+A1zns+s
A1x:
[2,2]
	0.000000	1000.000000
	0.000000	0.000000
A1zs:
[2,3]
	1.000000	-1.000000	-1000.000000
	500.000000	-500.000000	0.000000
A1zns:
A1s:
[2,1]
	1000.000000
	0.000000
print cap i coefs c1:
	0.000000	0.000000	0.000000	1000.000000	1.000000	-1.000000	-1000.000000	1000.000000
