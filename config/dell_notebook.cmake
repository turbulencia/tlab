# Dell Latitude - Notebook

if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE {PARALLEL})
endif()

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build
   set(ENV{FC} mpif90)
   set(CMAKE_Fortran_COMPILER mpif90)
   set(USER_Fortran_FLAGS "-cpp -std=legacy -ffree-form -ffree-line-length-none -fno-automatic")
   set(USER_Fortran_FLAGS_RELEASE "-ffpe-summary=none -O3 -fconvert=little-endian -O3 -ffast-math -ffinite-math-only -mtune=native -march=native") # -fallow-argument-mismatch
   add_definitions(-DUSE_FFTW -DUSE_MPI -DUSE_ALLTOALL -DUSE_MPI_IO)# -DIBM_DEBUG)
   set(CMAKE_BUILD_TYPE RELEASE)

else() # compiler for serial build
   set(ENV{FC} gfortran-10)
   set(CMAKE_Fortran_COMPILER gfortran-10)
   set(USER_Fortran_FLAGS "-cpp -std=legacy -ffree-form -ffree-line-length-none -fno-automatic -fallow-argument-mismatch")
   add_definitions(-DUSE_FFTW)# -DIBM_DEBUG)

   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
     set(USER_Fortran_FLAGS_RELEASE "-fconvert=big-endian -ffpe-summary=none -O3 -ffast-math -ffinite-math-only -mtune=native -march=native -funroll-loops")
     set(CMAKE_BUILD_TYPE RELEASE)

   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" )
     set(USER_Fortran_FLAGS_RELEASE "-fconvert=little-endian -ffpe-summary=none -O3 -ffast-math -ffinite-math-only -mtune=native -march=native -funroll-loops")
     set(CMAKE_BUILD_TYPE RELEASE)

   elseif( ${BUILD_TYPE} STREQUAL "DEBUG" )
     # set(USER_Fortran_FLAGS_DEBUG "-O0 -p -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal")
     set(USER_Fortran_FLAGS_DEBUG "-Og -ggdb3 -Wall -fbacktrace -fconvert=little-endian -ffpe-trap=invalid,zero,overflow")
     add_definitions(-D_DEBUG)
     set(CMAKE_BUILD_TYPE DEBUG)

   endif()

endif()

set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/usr/local/include")
set(FFTW_LIB           "/usr/lib/x86_64-linux-gnu/libfftw3.a")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

add_definitions(-DUSE_NETCDF)
set(NCDF_INCLUDE_DIR   "/usr/include")
set(NCDF_LIB           "-L/usr/lib/x86_64-linux-gnu -lnetcdff -lnetcdf")
set(INCLUDE_DIRS ${INCLUDE_DIRs} ${NCDF_INCLUDE_DIR})
set(LIBS ${LIBS} ${NCDF_LIB})