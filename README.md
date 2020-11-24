# Siconos/ACEF

Automatic circuit equation formulation 



**WORK IN PROGRESS**

Old simulation code (tag:acef_v1.0) works on ubuntu 10.04.

Start to update with a recent siconos


## Licence

### Siconos

Siconos is a program dedicated to modeling, simulation and control  of non smooth dynamical systems.

Copyright 2016 INRIA.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 
 See the License for the specific language governing permissions and limitations under the License.

### Ngspice

Ngspice is a mixed-level/mixed-signal circuit simulator. Its code
is based on three open source software packages: Spice3f5, Cider1b1
and Xspice.

The ngspice source code has evolved over time by integrating contributions
from various sources (e.g. Spice3f5, XSPICE, CIDER, numparam, tclspice and
others). Thus a mixture of license statements prevails.

The ngspice license is the `Modified BSD' license. This is adopted for all of
its source code, test and example files except for the files listed below.

## Directory

* docs/
* Examples/
* Sandbox/
* Tests/
* web/

* src/ This directory contains the code source of the automatic circuit
equation formulation (ACEF). It is based on Siconos.

  * `src/Spiceparser`: It builds the library libSiconosSpiceparser. It contains the netlist parser. 
  * `src/ACE`: It builds the library libSiconosACEF. It is linked with libSiconosSpiceparser. It contains the automatic circuit equation formulation.
  * `src/EXE`: It builds executables to run the simulation.

## Compilation

Depends on :

* Siconos (branch noselecs_v0) to be installed, see https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/index.html

* Spicelib from ngspice (copied in ./Spiceparser) see http://ngspice.sourceforge.net.  Do not identified the ngspice version (< 30).

* boost-graph, libXml2, blas/lapack

```
cd trunk/src/Spiceparser/
ln -s cmake SOURCE_SICONOS/Numerics/cmake
mkdir build
cd build
cmake ..
make install

cd trunk/src/ACE/
ln -s cmake SOURCE_SICONOS/Numerics/cmake
mkdir build
cd build
cmake ..
make install

cd trunk/src/EXE
ln -s cmake SOURCE_SICONOS/Numerics/cmake
mkdir build
cd build
cmake ..
make install
```

### Test
cd Tests/ACE
../../trunk/src/EXE/build/src/noselect LC.cir ENUM

