# JUWELS (FZ JUELICH)

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

#set(USER_Fortran_FLAGS          " -fpp ${USER_profile_FLAGS} -nbs -save-temps -heap-arrays -vec-threshold50  ${USER_omp_FLAGS} " )
#set(USER_Fortran_FLAGS_RELEASE  " -march=skylake-avx512 -axcommon-avx512,SSE4.2 -O3 -ipo" )
#set(USER_Fortran_FLAGS_RELEASE  " -axCORE-AVX2 -O3 -ipo" )

set(USER_Fortran_FLAGS         "-cpp -ffree-form -ffree-line-length-2048 -fno-automatic -fallow-argument-mismatch " )
set(USER_Fortran_FLAGS_RELEASE "-march=skylake-avx512 -mtune=skylake-avx512 -O3 -ffinite-math-only -fprefetch-loop-arrays --param prefetch-latency=300 " )

set(USER_Fortran_FLAGS_DEBUG    " -g -traceback -debug all " )

if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

add_definitions(-DUSE_FFTW -DUSE_BLAS -DUSE_MKL -DUSE_NETCDF)

#set(FFTW_INCLUDE_DIR "/gpfs/software/juwels/stages/2018a/software/FFTW/3.3.7-ipsmpi-2018a/include/")
#set(FFTW_LIB         "/gpfs/software/juwels/stages/2018a/software/FFTW/3.3.7-ipsmpi-2018a/lib/libfftw3.a")
set(NC_INCLUDE_DIR   "/p/software/juwels/stages/2024/software/netCDF-Fortran/4.6.1-gpsmpi-2023a/include/ " )
set(FFTW_LIB         "-lfftw3")
#set(BLAS_LIB         "/gpfs/software/juwels/stages/2018a/software/imkl/2018.2.199-iimpi-2018a/mkl/lib/intel64/libmkl_blas95_ilp64.a")
set(BLAS_LIB         "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(NCDF_LIB         "-lnetcdff")

set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NC_INCLUDE_DIR} )
set(LIBS         ${NCDF_LIB} ${FFTW_LIB} ${BLAS_LIB})

set(NC_LIB             "-L/usr/lib -lnetcdff -lnetcdf")
set(INCLUDE_DIRS ${INCLUDE_DIRS} )
set(LIBS ${LIBS} ${NC_LIB})

if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  add_definitions(-DUSE_PSFFT)

  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qopenmp")
  set(PSFFT_COMPILE_FLAGS "-qopenmp")

  set(INCLUDE_DIRS ${INCLUDE_DIRS} ${PSFFT_INCLUDE_DIR})

endif()

set(GNU_SED "gsed")
