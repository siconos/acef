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
	________U2_3	________I0_3
Zs
	________V0	________V1	________V2	________V3	________I1_0
Zns

---------------------------------------
equation	________U2_3	________I0_3	________U2_3	________I0_3	________V0	________V1	________V2	________V3	________I1_0
KCL2*0		-0.020000000000	0.000000000000	0.000000000000	0.000000000000	0.000000000000	0.010000000000	-0.010000000000	0.000000000000	0.000000000000	0.000000000000
IND*1		0.000000000000	0.030000000000	0.000000000000	0.000000000000	1.000000000000	0.000000000000	0.000000000000	-1.000000000000	0.000000000000	0.000000000000
KCL0		0.000000000000	0.000000000000	0.000000000000	-1.000000000000	0.000000000000	0.000000000000	0.000000000000	0.000000000000	1.000000000000	0.000000000000
VD		0.000000000000	0.000000000000	0.000000000000	0.000000000000	1.000000000000	-1.000000000000	0.000000000000	0.000000000000	0.000000000000	-3.000000000000
TEN		0.000000000000	0.000000000000	1.000000000000	0.000000000000	0.000000000000	0.000000000000	-1.000000000000	1.000000000000	0.000000000000	0.000000000000
A:
[2,2]
	-0.020000	0.000000
	0.000000	0.030000
inv A:
[2,2]
	-50.000000	0.000000
	0.000000	33.333333
B:
[2,2]
	0.000000	0.000000
	0.000000	0.000000
C:
[2,5]
	0.000000	0.010000	-0.010000	0.000000	0.000000
	1.000000	0.000000	0.000000	-1.000000	0.000000
D:
s:
[2,1]
	0.000000
	0.000000
system x'=A1x+A1zs+A1zns+s
A1x:
[2,2]
	0.000000	0.000000
	0.000000	0.000000
A1zs:
[2,5]
	0.000000	-0.500000	0.500000	0.000000	0.000000
	33.333333	0.000000	0.000000	-33.333333	0.000000
A1zns:
A1s:
[2,1]
	0.000000
	0.000000
print cap i coefs c1:
	0.000000	0.000000	0.000000	0.000000	0.000000	-0.500000	0.500000	0.000000	0.000000	0.000000
