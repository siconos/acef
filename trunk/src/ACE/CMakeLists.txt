cmake_minimum_required(VERSION 2.4)
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran) # see lapack macros

INCLUDE(SiconosProject)
INCLUDE(FindPkgConfig)

COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)
SET(Boost_FIND_COMPONENTS graph)
COMPILE_WITH(Boost REQUIRED)
IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.5)
   if (NOT BOOST_GRAPH_FOUND)
      message( "boost graph must be installed")
   endif (NOT BOOST_GRAPH_FOUND)
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.5)

COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(LAPACK REQUIRED)
COMPILE_WITH(BLAS REQUIRED)
COMPILE_WITH(ATLAS)
IF (NOT ATLAS_INCLUDE_PATH)
 SET(ATLAS_INCLUDE_PATH ${SiconosKernel_INCLUDE_DIRS}/AtlasLocal)
 include_directories(${ATLAS_INCLUDE_PATH})

ENDIF (NOT ATLAS_INCLUDE_PATH)
SET(PACKAGE_DESCRIPTION_SUMMARY 
  "ACEF is a Automatic Circuit Equation Formution library.")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
ACEF is a Automatic Circuit Equation Formution library.")

SICONOS_PROJECT(ACEF 3 0 0)
