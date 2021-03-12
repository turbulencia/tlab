cmake_minimum_required (VERSION 2.8...3.13)
# specify the location of loadable modules (.cmake extension)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../config)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# load system specific settings if not set, force default.cmake
# cmake ../src -DSYST=mpipc
if(NOT SYST)
  set(SYST default)
endif()
include(${SYST} OPTIONAL RESULT_VARIABLE SYSTINC)

# trigger fatal error if illegal module is loaded
if(${SYSTINC} STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Config file config/" ${SYST} ".cmake does not exist.")
endif()

# set the default build type to RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
endif()

# only set the compiler flags when the cache is created
# to enable editing of the flags in the CMakeCache.txt file
if(NOT HASCACHE)
  set(CMAKE_Fortran_FLAGS ${USER_Fortran_FLAGS} CACHE STRING
    "Flags used by the compiler during all build types." FORCE)
  set(CMAKE_Fortran_FLAGS_DEBUG ${USER_Fortran_FLAGS_DEBUG} CACHE STRING
    "Flags used by the compiler during debug builds." FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} CACHE STRING
    "Flags used by the compiler during release builds." FORCE)
  message(STATUS "Build Type: " ${CMAKE_BUILD_TYPE})
  set(HASCACHE TRUE CACHE BOOL "CMakeCache.txt created." FORCE)
  # make sure that ccmake only contains build type and MPI setting
  mark_as_advanced(HASCACHE)
  mark_as_advanced(CMAKE_INSTALL_PREFIX)
endif()

# print the compiler flags to the screen
if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  message(STATUS "Compiler flags: " ${CMAKE_Fortran_FLAGS} " " ${CMAKE_Fortran_FLAGS_RELEASE})
else()
  message(STATUS "Compiler flags: " ${CMAKE_Fortran_FLAGS} " " ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

# start the project only after the system specific settings are loaded
project(tlab Fortran)

# set the include dirs
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${PSFFT_INCLUDE_DIR})

# produce mod files
add_subdirectory(modules)
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/modules)

# produce libraries
add_subdirectory(utils)
add_subdirectory(fdm)
add_subdirectory(linear)
add_subdirectory(filters)
add_subdirectory(averages)
add_subdirectory(pdfs)
add_subdirectory(thermodynamics)
add_subdirectory(fitpack)
add_subdirectory(quadpack)
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
add_subdirectory(mpi)
endif()
add_subdirectory(operators)
add_subdirectory(io)
add_subdirectory(mappings)
add_subdirectory(lagrange)

# produce executables
add_subdirectory(tools)

# produce validation tools
add_subdirectory(valid)

#add_custom_target(examples message(STATUS "We are making examples"))
