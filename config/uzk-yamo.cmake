# University of Colgne - Workstation (yamo) [Ubuntu 22.04 LTS. 64-bit]
if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE {PRALLEL})
endif()
if ( NOT HYBRID )
   set(HYBRID FALSE)
else()
   message(WARNING "Compiling for hybrid openMP/MPI usage")
endif()
if ( NOT PROFILE )
   set(PROFILE FALSE)
endif()
if ( ${PROFILE} STREQUAL "TRUE" )
   set(USER_profile_FLAGS "-g")
endif()

# compiler for parallel build
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" )
   set(ENV{FC} mpif90)
   set(CMAKE_Fortran_COMPILER mpif90)
   set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic -fallow-argument-mismatch")
   set(USER_Fortran_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")
   add_definitions(-DUSE_FFTW -DUSE_MPI -DUSE_MPI_IO)
   set(CMAKE_BUILD_TYPE RELEASE)

# compiler for serial build
elseif( ${BUILD_TYPE} STREQUAL "SERIAL" )
   set(ENV{FC} gfortran-11)
   set(CMAKE_Fortran_COMPILER gfortran)
   set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic")
   add_definitions(-DUSE_FFTW)
   set(USER_Fortran_FLAGS_RELEASE "-ffpe-summary=none -O3 -ffast-math -mtune=native -march=native -fallow-argument-mismatch")
   set(CMAKE_BUILD_TYPE RELEASE)

# compiler for debug build
elseif( ${BUILD_TYPE} STREQUAL "DEBUG" )
   set(ENV{FC} gfortran-11)
#  set(USER_Fortran_FLAGS_DEBUG "-O0 -p -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal")
   set(USER_Fortran_FLAGS_DEBUG "-cpp -ffree-line-length-none -O0 -ggdb -Wall -fbacktrace -fconvert=big-endian -ffpe-trap=invalid,zero,overflow -fallow-argument-mismatch")
   add_definitions(-D_DEBUG -DUSE_FFTW -DIBM_DEBUG -DTRACE_ON)
   set(CMAKE_BUILD_TYPE DEBUG)
endif()

set(GNU_SED "gsed")

# NETCDF (HAWK for dns_transpose block size)
add_definitions(-DUSE_NETCDF -DHLRS_HAWK) 
set(NC_LIB "-L/usr/lib -lnetcdff -lnetcdf")
# FFTW
set(FFTW_INCLUDE_DIR   "/usr/include/")
set(FFTW_LIB           "/usr/lib/x86_64-linux-gnu/libfftw3.a")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB} ${NC_LIB})
