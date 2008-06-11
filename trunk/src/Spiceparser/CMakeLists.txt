SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

ENABLE_LANGUAGE(C)

INCLUDE(SiconosProject)

SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Spiceparser is a netlist parser library.")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/Spiceparser is a netlist parser library.")

SICONOS_PROJECT(Spiceparser 3 0 0)
