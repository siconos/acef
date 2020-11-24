#
# Here comes Spiceparser
#
PROJECT(SiconosSpiceparser C)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosSpiceparser_VERSION "${VERSION}")

#
# Published headers
#
# if <PROJECT>_HDRS unset -> all headers under <PROJECT>_DIRS

SET(SiconosSpiceparser_HDRS
	Spiceparser.h)
	
#
# Project directories 
#
SET(SiconosSpiceparser_DIRS
  ..
  .
  frontend
  frontend/help
  frontend/numparam
  frontend/parser
  frontend/plotting
  frontend
  frontend/wdisp
  include
  maths
  maths/cmaths
  maths/deriv
  maths/misc
  maths/ni
  maths/poly
  maths/sparse
  misc
  spicelib
  spicelib/analysis
  spicelib/devices
  spicelib/devices/asrc
  spicelib/devices/bjt
  spicelib/devices/bjt2
  spicelib/devices/bsim1
  spicelib/devices/bsim2
  spicelib/devices/bsim3
  spicelib/devices/bsim3soi
  spicelib/devices/bsim3soi_dd
  spicelib/devices/bsim3soi_fd
  spicelib/devices/bsim3soi_pd
  spicelib/devices/bsim3v0
  spicelib/devices/bsim3v1
  spicelib/devices/bsim3v1a
  spicelib/devices/bsim3v1s
  spicelib/devices/bsim4
  spicelib/devices/cap
  spicelib/devices/cccs
  spicelib/devices/ccvs
  spicelib/devices/comp
  spicelib/devices/cpl
  spicelib/devices/csw
  spicelib/devices/dio
  spicelib/devices/hfet1
  spicelib/devices/hfet2
  spicelib/devices/hisim
  spicelib/devices/ind
  spicelib/devices/isrc
  spicelib/devices/jfet
  spicelib/devices/jfet2
  spicelib/devices/ltra
  spicelib/devices/mes
  spicelib/devices/mesa
  spicelib/devices/mos1
  spicelib/devices/mos2
  spicelib/devices/mos3
  spicelib/devices/mos6
  spicelib/devices/mos9
  spicelib/devices/nbjt
  spicelib/devices/nbjt2
  spicelib/devices/numd
  spicelib/devices/numd2
  spicelib/devices/res
  spicelib/devices/soi3
  spicelib/devices/sw
  spicelib/devices/tra
  spicelib/devices/txl
  spicelib/devices/urc
  spicelib/devices/vbic
  spicelib/devices/vccs
  spicelib/devices/vcvs
  spicelib/devices/vsrc
  spicelib/parser
  )


APPEND_FLAGS(CMAKE_C_FLAGS "-DSIMULATOR ")
#
# Installation layout
#
SET(SiconosSpiceparser_INSTALL_INCLUDE_DIR include/Siconos/Spiceparser)

# SpiceparserConfig.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/SpiceparserConfig.h DESTINATION include/Siconos/Spiceparser)

SET(SiconosSpiceparser_LINKER_LANGUAGE C)

LIBRARY_PROJECT_SETUP()

#TARGET_LINK_LIBRARIES(${PROJECT_NAME}_static "-lm -lXaw -lXmu")
#TARGET_LINK_LIBRARIES(${PROJECT_NAME}_shared "-lm -lXaw -lXmu")


IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)