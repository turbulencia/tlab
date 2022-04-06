# MacBook

if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE PARALLEL)
endif()


if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build
   set(ENV{FC} /opt/local/bin/mpif90-mpich-clang11)
   set(CMAKE_Fortran_COMPILER /opt/local/bin/mpif90-mpich-clang11)
   set(USER_Fortran_FLAGS     "-fallow-argument-mismatch -ffree-form -ffree-line-length-2048 -fno-automatic -O3 -cpp -arch x86_64 -ffast-math -ffinite-math-only -funroll-loops -mtune=native -Wunused")
#   set(USER_Fortran_FLAGS_RELEASE  "-O3 -cpp -arch x86_64 -ffast-math -ffinite-math-only -funroll-loops -mtune=native ")

   add_definitions(-DUSE_MPI -DUSE_ALLTOALL -DUSE_MPI_IO -DUSE_FFTW)

   set(CMAKE_BUILD_TYPE PARALLEL )

else() # compiler for serial build
   set(ENV{FC} /opt/local/bin/gfortran-mp-10)
   set(USER_Fortran_FLAGS          "-cpp -fallow-argument-mismatch -ffree-form -ffree-line-length-2048 -fno-automatic")

#   add_definitions(-DUSE_PSFFT -DUSE_FFTW)
   add_definitions( -DUSE_FFTW)

   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -fconvert=big-endian    -mtune=native -ffast-math -ffinite-math-only -funroll-loops")
     set(CMAKE_BUILD_TYPE RELEASE)

   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" )
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -fconvert=little-endian -mtune=native -ffast-math -ffinite-math-only -funroll-loops")
#     add_definitions(-DTRACE_ON)

   else()
     set(USER_Fortran_FLAGS_DEBUG    "-Og -ggdb3 -Wall -Wno-unknown-pragmas")
     add_definitions(-D_DEBUG)
     set(CMAKE_BUILD_TYPE DEBUG)
   endif()


endif()

add_definitions(-DUSE_NETCDF)

set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/opt/local/include")
set(FFTW_LIB           "/opt/local/lib/libfftw3.a")
set(NCDF_INCLUDE_DIR   "/opt/local/include")
set(NCDF_LIBPATH           "/opt/local/lib")
#set(FFTW_INCLUDE_DIR   "/usr/local/fftw_intel14/include")
#set(FFTW_LIB           "/usr/local/fftw_intel14/lib/libfftw3.a")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NCDF_INCLUDE_DIR} )
set(LIBS ${FFTW_LIB} ${NCDF_LIBPATH}/libnetcdff.a ${NCDF_LIBPATH}/libnetcdf.dylib )

add_definitions(-DRESTRICTKEYWORD=__restrict__)


if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qsmp=omp")
  set(PSFFT_COMPILE_FLAGS "-qsmp=omp:noauto:nested_par")

  set(INCLUDE_DIRS ${INCLUDE_DIRS} ${PSFFT_INCLUDE_DIR})

  add_definitions(-WF,-DUSE_PSFFT)
endif()
