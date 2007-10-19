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
x
	________U2_0
Zs
	________V0	________V1	________V2	________I0_1
Zns

---------------------------------------
equation	________U2_0	________U2_0	________V0	________V1	________V2	________I0_1
KCL0		0.000000000000	0.000000000000	0.000000000000	0.000000000000	0.000000000000	-1.000000000000	0.000000000000
KCL1		0.000000000000	0.000000000000	0.000000000000	-0.010000000000	0.010000000000	1.000000000000	0.000000000000
KCL2*		-0.020000000000	0.000000000000	0.000000000000	0.010000000000	-0.010000000000	0.000000000000	0.000000000000
VD		0.000000000000	0.000000000000	-1.000000000000	1.000000000000	0.000000000000	0.000000000000	-3.000000000000
TEN		0.000000000000	1.000000000000	1.000000000000	0.000000000000	-1.000000000000	0.000000000000	0.000000000000
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
