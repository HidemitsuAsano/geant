#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(knucl)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# include root
if( APPLE )
SET(CMAKE_MODULE_PATH /usr/local/geant4/geant4.9.5.p01/cmake/Modules/)
else( APPLE )
if( $ENV{HOSTNAME} STREQUAL Ag.riken.jp )
SET(CMAKE_MODULE_PATH /w/e15/common/GEANT/geant4.10.00_source/cmake/Modules/)
else( $ENV{HOSTNAME} STREQUAL Ag.riken.jp )
SET(CMAKE_MODULE_PATH /group/had/knucl/knucl_soft/geant4.9.5.p02/cmake/Modules/)
endif( $ENV{HOSTNAME} STREQUAL Ag.riken.jp )
endif( APPLE )

find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIR})
find_package(ROOT REQUIRED)
#include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# set root dictionary
#set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lEve)
#target_link_libraries(knucl ${ROOT_LIBRARIES} Eve)
#set(CMAKE_EXE_LINKER_FLAGS "-L/usr/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic")
#ROOT_GENERATE_DICTIONARY(EventDict ${CMAKE_CURRENT_SOURCE_DIR}/Event.h LINKDEF EventLinkDef.h)
#ROOT_LINKER_LIBRARY(Event Event.cxx EventDict.cxx LIBRARIES Hist MathCore)
#ROOT_EXECUTABLE(eventexe MainEvent.cxx LIBRARIES Event RIO Tree Hist Net)
#ROOT_ADD_TEST(test-event COMMAND eventexe)

#ROOT_GENERATE_DICTIONARY(KnuclRootDataCint.cc ${PROJECT_SOURCE_DIR}/include/KnuclRootData.h LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)
#ROOT_LINKER_LIBRARY(KnuclRootData ${PROJECT_SOURCE_DIR}/src/KnuclRootData.cc KnuclRootDataCint.cc LIBRARIES Hist MathCore)

EXECUTE_PROCESS(COMMAND rootcint -f KnuclRootDataCint.cc -c ${PROJECT_SOURCE_DIR}/include/KnuclRootData.h ${PROJECT_SOURCE_DIR}/include/ComCrossSectionTable.hh ${PROJECT_SOURCE_DIR}/include/LinkDef.h)
EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_EXE_LINKER_FLAGS "${ROOT_LD_FLAGS} KnuclRootDataCint.cc")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(knucl knucl.cc ${sources} ${headers})
target_link_libraries(knucl ${Geant4_LIBRARIES})
target_link_libraries(knucl ${ROOT_LIBRARIES})
#target_link_libraries(knucl KnuclRootData)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
	.rootlogon.C
	CS.list
	K3He.C
	KnuclSetting.card
	KppShape.root
	Magnet_Field_Doraemon_3.txt
	dump.mac
	nophys.mac
	osc.dat
	pdg_table.txt
	run.mac
	ushiwaka_field_map.txt
	vis.mac
	vis_qt.mac
	vrml.mac
  )
  
foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(knucl3_ag DEPENDS knucl)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS knucl DESTINATION bin)

EXECUTE_PROCESS(COMMAND ln -s ${PROJECT_SOURCE_DIR}/src/KnuclRootData.cc .)
EXECUTE_PROCESS(COMMAND ln -s ${PROJECT_SOURCE_DIR}/include/KnuclRootData.h .)
EXECUTE_PROCESS(COMMAND ln -s ${PROJECT_SOURCE_DIR}/src/ComCrossSectionTable.cc .)
EXECUTE_PROCESS(COMMAND ln -s ${PROJECT_SOURCE_DIR}/include/ComCrossSectionTable.hh .)
EXECUTE_PROCESS(COMMAND ln -s ${PROJECT_SOURCE_DIR}/Makefile.lib .)
