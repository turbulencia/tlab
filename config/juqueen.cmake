########################################################
#  IBM BlueGene/Q (JUQUEEN)
#
# If SMP, use -qsmp=omp at compiling and linking, 
# remove the -qsave in FFLAGS, and add -DUSE_OPENMP to
# the CPPFLAGS
#
# Add -DUSE_BLAS to CPPFLAGS if you want to use the BLAS
# routines from the ESSL library
#
# Add -DUSE_ESSL to CPPFLAGS if you want to use the
# transpose routine from the ESSL library
#
########################################################
if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE PARALLEL)
endif()
 
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # compiler for parallel build
   set(ENV{FC} /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf95_r)
#   set(CMAKE_Fortran_COMPILER mpixlf95_r) 
   set(USER_Fortran_FLAGS         "-qnoescape ")
   set(USER_Fortran_FLAGS_RELEASE "-O3 -qmaxmem=-1 -qarch=qp -qtune=qp -qxflag=ngenstub -qhot -qalias=noaryovrlp:nopteovrlp -qassert=contiguous ")
   set(USER_Fortran_FLAGS_DEBUG   "-O0 -pg -qfullpath")

   add_definitions(-d -qsuffix=cpp=f90 -WF,-DUSE_MPI,-DUSE_MPI_IO,-DUSE_FFTW,-DUSE_ESSL)

   set(CMAKE_BUILD_TYPE RELEASE)

else() # compiler for serial build; not developed
  message( FATAL_ERROR "Only PARALLEL mode developed. CMake will exit." )
endif()

set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/bgsys/local/fftw3/3.3.2/fftw_g/include")
set(FFTW_LIB           "/bgsys/local/fftw3/3.3.2/fftw_g/lib/libfftw3.a")
set(FFTWTHREADS_LIB    "/bgsys/local/fftw3/3.3.2/fftw_g/lib/libfftw3_threads.a")
set(ESSL_LIB           "/bgsys/local/lib/libesslbg.a")

set(LIBS ${FFTW_LIB} ${ESSL_LIB})

if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qsmp=omp")
  set(PSFFT_COMPILE_FLAGS "-qsmp=omp:noauto:nested_par")
  
  add_definitions(-WF,-DUSE_PSFFT)
endif()
