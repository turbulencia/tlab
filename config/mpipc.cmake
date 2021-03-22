if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE BIG)
endif()

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build
   set(ENV{FC} mpif90)
   set(CMAKE_Fortran_COMPILER mpif90)
   set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic")
   set(USER_Fortran_FLAGS_RELEASE "-fconvert=little-endian -O3 -ffast-math -mtune=native -march=native")

   add_definitions(-DUSE_FFTW -DUSE_MPI -DUSE_MPI_IO)

   set(CMAKE_BUILD_TYPE RELEASE)

else() # compiler for serial build
   set(ENV{FC} gfortran)
   set(CMAKE_Fortran_COMPILER gfortran)
   set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic")

   add_definitions(-DUSE_FFTW)

   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
     set(USER_Fortran_FLAGS_RELEASE "-fconvert=big-endian -ffpe-summary=none -O3 -ffast-math -mtune=native -march=native")

     set(CMAKE_BUILD_TYPE RELEASE)

   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" )
     set(USER_Fortran_FLAGS_RELEASE "-fconvert=little-endian -ffpe-summary=none -O3 -ffast-math -mtune=native -march=native")

     set(CMAKE_BUILD_TYPE RELEASE)

   else()
#     set(USER_Fortran_FLAGS_DEBUG "-O0 -p -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal")
     set(USER_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -fbacktrace -fconvert=big-endian -ffpe-trap=invalid,zero,overflow")

     add_definitions(-D_DEBUG)

     set(CMAKE_BUILD_TYPE DEBUG)

   endif()

endif()

set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/sw/jessie-x64/fftw-3.3.4-openmp-gccsys/include")
set(FFTW_LIB           "/sw/jessie-x64/fftw-3.3.4-openmp-gccsys/lib/libfftw3.a")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB} )
