# JURECA (FZ JUELICH)

if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE PARALLEL)
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

# compiler for parallel build and hybrid flags
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
   set(ENV{FC} scorep-mpif90)
   if ( ${PROFILE} STREQUAL "TRUE" )
      set(CMAKE_Fortran_COMPILER scorep-mpif90)
   else ()
      set(CMAKE_Fortran_COMPILER mpif90)
   endif()

   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL)

   if ( ${HYBRID} STREQUAL "TRUE" )
      set(USER_omp_FLAGS " -qopenmp " )
   endif()

# compiler for serial build
else( ${BUILD_TYPE} STREQUAL "SERIAL" )
   set(ENV{FC} ifort)
   set(CMAKE_Fortran_COMPILER ifort)
endif()

set(USER_Fortran_FLAGS          " -fpp ${USER_profile_FLAGS} -O3 -nbs -save-temps -heap-arrays -simd -vec-threshold50 -unroll-aggressive ${USER_omp_FLAGS} " )
set(USER_Fortran_FLAGS_RELEASE  " -march=native -qopt-prefetch -O3 -ipo" )
set(USER_Fortran_FLAGS_DEBUG    " -g -traceback -debug all ")

if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

add_definitions(-DUSE_FFTW -DUSE_BLAS)

set(FFTW_LIB         "-lfftw3")
set(BLAS_LIB         "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS             "${FFTW_LIB} ${BLAS_LIB}"  )

set(GNU_SED "gsed")

if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  add_definitions(-DUSE_PSFFT)

  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qopenmp")
  set(PSFFT_COMPILE_FLAGS "-qopenmp")

  set(INCLUDE_DIRS ${INCLUDE_DIRS} ${PSFFT_INCLUDE_DIR})

endif()
