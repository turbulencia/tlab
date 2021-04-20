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
   set(ENV{FC} mpif90) 

   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL)

   if ( ${HYBRID} STREQUAL "TRUE" )
     set(USER_omp_FLAGS " -fopenmp " )
     add_definitions(-DUSE_OPENMP) 
   endif()

# compiler for serial build
else( ${BUILD_TYPE} STREQUAL "SERIAL" ) 
  set(ENV{FC} mpif90)
  set(CMAKE_Fortran_COMPILER mpif90)

endif()     


set(DRAGONEGG_FLAGS "-finline-aggressive -fslp-vectorize  -fmerge-all-constants") #  -mmadd4 -mfp64 -enable-strided-vectorization")

set(USER_Fortran_FLAGS         "-cpp ${USER_omp_FLAGS}") #-fallow-argument-mismatch from gnu-version10
set(USER_Fortran_FLAGS_RELEASE "-O3 -march=znver2 -mtune=znver2 ${DRAGONEGG_FLAGS}") #these will be ignored:  -fprefetch-loop-arrays --param prefetch-latency=300") 
set(USER_Fortran_FLAGS_DEBUG   "-O0 -g -debug -ffpe-trap=all") 

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

add_definitions(-DUSE_FFTW -DHLRS_HAWK) # -DUSE_BLAS -DUSE_MKL) 

set(FFTW_LIB   "-lfftw3")
#set(FFTW_INCLUDE_DIR "/opt/hlrs/spack/rev-004_2020-06-17/fftw/3.3.8-clang-9.0.0-2buapgdw/include/")
set(LIBS             "${FFTW_LIB} -l:libamdlibm.a -lm")

set(GNU_SED "gsed")

if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  add_definitions(-DUSE_PSFFT)

  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qopenmp")
  set(PSFFT_COMPILE_FLAGS "-qopenmp")

endif()
