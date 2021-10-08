if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE BIG)
endif()

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build
   set(ENV{FC} mpif90)
   set(CMAKE_Fortran_COMPILER mpif90)
   set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic")
   set(USER_Fortran_FLAGS_RELEASE "-O3 -ffast-math -mtune=native -march=native")

   add_definitions(-DUSE_FFTW -DUSE_MPI -DUSE_MPI_IO)

   set(CMAKE_BUILD_TYPE RELEASE)

elseif( ${BUILD_TYPE} STREQUAL "SERIAL" ) # compiler for serial build
   set(ENV{FC} gfortran)
   set(CMAKE_Fortran_COMPILER gfortran)
   set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic")

   add_definitions(-DUSE_FFTW)

   set(USER_Fortran_FLAGS_RELEASE "-ffpe-summary=none -O3 -ffast-math -mtune=native -march=native")

   set(CMAKE_BUILD_TYPE RELEASE)

elseif( ${BUILD_TYPE} STREQUAL "DEBUG" )
#     set(USER_Fortran_FLAGS_DEBUG "-O0 -p -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal")
     set(USER_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -fbacktrace -fconvert=big-endian -ffpe-trap=invalid,zero,overflow")

     add_definitions(-D_DEBUG)

     set(CMAKE_BUILD_TYPE DEBUG)

endif()

add_definitions(-DUSE_NETCDF -DHLRS_HAWK) 
set(NC_LIB "-L/usr/lib -lnetcdff -lnetcdf")


set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/usr/include/")
set(FFTW_LIB           "/usr/lib/x86_64-linux-gnu/libfftw3.a")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB} ${NC_LIB})