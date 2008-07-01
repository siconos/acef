#
# Here comes ACEF_EXE
#
PROJECT(SiconosACEF_EXE CXX)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosACEF_EXE_VERSION "${VERSION}")

#
# Published headers
#
# if <PROJECT>_HDRS unset -> all headers under <PROJECT>_DIRS

SET(SiconosACEF_EXE_HDRS
	ACEF.h)
	
#
# Project directories 
#
SET(SiconosACEF_EXE_DIRS
  ..
  .)


#
# Installation layout
#
SET(SiconosACEF_EXE_INSTALL_INCLUDE_DIR include/Siconos/ACEF_EXE)

# ACEF_EXEConfig.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/ACEF_EXEConfig.h DESTINATION include/Siconos/ACEF_EXE)

SET(SiconosACEF_EXE_LINKER_LANGUAGE CXX)
add_executable (acef main.cpp)
add_executable (noselect mainSiconos.cpp)

target_link_libraries (acef SiconosACEF SiconosNumerics SiconosKernel SiconosSpiceparser)
target_link_libraries (noselect SiconosACEF SiconosNumerics SiconosKernel SiconosSpiceparser)

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)