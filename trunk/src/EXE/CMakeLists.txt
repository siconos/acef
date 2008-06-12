SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran) # see lapack macros

INCLUDE(SiconosProject)
COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)
COMPILE_WITH(SiconosACEF REQUIRED)
COMPILE_WITH(Boost REQUIRED)
COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(LAPACK REQUIRED)
COMPILE_WITH(BLAS REQUIRED)
COMPILE_WITH(ATLAS)


SET(PACKAGE_DESCRIPTION_SUMMARY 
  "ACEF_EXE is a Automatic Circuit Equation Formution software based on ACEF.")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
ACEF is a Automatic Circuit Equation Formution software based on ACEF.")

SICONOS_PROJECT(ACEF_EXE 3 0 0)
