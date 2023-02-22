#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
cmake_host_system_information(RESULT host QUERY HOSTNAME)
cmake_host_system_information(RESULT ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
if(host MATCHES ".kek\.jp")
  message(STATUS "The hostname is ${host}")
  set(CMAKE_C_COMPILER /opt/gcc-8.3.0/bin/gcc)
  set(CMAKE_CXX_COMPILER /opt/gcc-8.3.0/bin/g++)
endif()
project(HYPTPC)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build HypTPC with Geant4 UI and Vis drivers" ON)
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

#-- Additional
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

#--

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(g4hyptpc main.cc ${sources} ${headers})
target_link_libraries(g4hyptpc ${Geant4_LIBRARIES})
target_link_libraries(g4hyptpc ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build HYPTPC. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(HYPTPC_SCRIPTS
  test.mac
  vis.mac
  )

foreach(_script ${HYPTPC_SCRIPTS})
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
add_custom_target(HYPTPC DEPENDS g4hyptpc)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS g4hyptpc DESTINATION bin)
