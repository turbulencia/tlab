# SUPERMUC

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
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" )
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

#set(USER_Fortran_FLAGS          " -fpp ${USER_profile_FLAGS} -nbs -save-temps -simd-vec-threshold50 -unroll-aggressive ${USER_omp_FLAGS} " )
set(USER_Fortran_FLAGS          " -fpp ${USER_profile_FLAGS} -O3 -nbs -save-temps -heap-arrays -simd -vec-threshold50 -unroll-aggressive ${USER_omp_FLAGS} " )
set(USER_Fortran_FLAGS_RELEASE  " -march=skylake-avx512 -axcommon-avx512,SSE4.2 -qopt-prefetch -O3 -ipo" )
#set(USER_Fortran_FLAGS_RELEASE  " -axCORE-AVX2 -qopt-prefetch -O3 -ipo" )
set(USER_Fortran_FLAGS_DEBUG    " -g -traceback -debug all ")

if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

add_definitions(-DUSE_BLAS -DUSE_MKL)
set(BLAS_LIB         "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(LIBS         ${BLAS_LIB})

add_definitions(-DUSE_FFTW)
set(FFTW_INCLUDE_DIR "/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/skylake_avx512/fftw/3.3.8-intel-ysseeml/include")
#set(FFTW_LIB         "/gpfs/software/juwels/stages/2018a/software/FFTW/3.3.7-ipsmpi-2018a/lib/libfftw3.a")
set(FFTW_LIB         $ENV{FFTW_LIB})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${FFTW_INCLUDE_DIR})
set(LIBS         ${LIBS}         ${FFTW_LIB})

add_definitions(-DUSE_NETCDF)
set(NC_INCLUDE_DIR     "/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/skylake_avx512/netcdf-fortran/4.5.3-intel-k5dpfai/include" "/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/skylake_avx512/netcdf-c/4.7.4-intel-zhencv3/include")
set(NC_LIB             $ENV{NETCDF_FORTRAN_SHLIB})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${NC_INCLUDE_DIR})
set(LIBS         ${LIBS}         ${NC_LIB})

set(GNU_SED "gsed")
