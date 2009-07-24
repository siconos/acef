#
# Here comes ACEF
#
PROJECT(SiconosACEF CXX)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosACEF_VERSION "${VERSION}")

#
# Published headers
#
# if <PROJECT>_HDRS unset -> all headers under <PROJECT>_DIRS

SET(SiconosACEF_HDRS
ACEF.h
ace.h
aceMatrix.h
acetime.h
aceVector.h
algo.h
componentarb.h
componentbjt.h
componentcap.h
componentcomp.h
componentdio.h
componentdyn.h
component.h
componentind.h
componentisrc.h
componentlinear.h
componentmos.h
componentnlinear.h
componentres.h
componentvccs.h
componentvcvs.h
componentvsrc.h
equationcap.h
equation.h
equation_nl.h
equationind.h
equationkcl.h
equationten.h
equationvd.h
extern.h
graph.h
linearsystem.h
linearsystemDAE.h
mlcp.h
unknown.h
acefRelation.h)

#
# Project directories 
#
SET(SiconosACEF_DIRS
  ..
  .)


#
# Installation layout
#
SET(SiconosACEF_INSTALL_INCLUDE_DIR include/Siconos/ACEF)

# ACEFConfig.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/ACEFConfig.h DESTINATION include/Siconos/ACEF)

SET(SiconosACEF_LINKER_LANGUAGE CXX)

LIBRARY_PROJECT_SETUP()

#TARGET_LINK_LIBRARIES(${PROJECT_NAME}_static "-lm -lXaw -lXmu")
#TARGET_LINK_LIBRARIES(${PROJECT_NAME}_shared "-lm -lXaw -lXmu")


IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)