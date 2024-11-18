# HLRS HUNTER-GPU System using new AMD compiler (Stuttgart / test at AMD- RAC-PLANO) 

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

if ( NOT ACCELERATE )
  set(ACCELERATE "FALSE")
elseif( ${ACCELERATE} STREQUAL "TRUE" )
  set(USER_accel_FLAGS "-fopenmp --offload-arch=gfx942 -cpp")
  add_definitions(-DUSE_APU)
endif() 

if ( ${PROFILE} STREQUAL "TRUE" )  
   set(USER_profile_FLAGS "-g")
endif() 

# compiler for parallel build and hybrid flags	  
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
   set(ENV{FC} /shared/apps/rhel8/opt/rocm-afar-5891/bin/amdflang-new ) 

   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL)

   if ( ${HYBRID} STREQUAL "TRUE" )
     set(USER_omp_FLAGS " -fopenmp " )
     add_definitions(-DUSE_OPENMP) 
   endif()

# compiler for serial build
else( ${BUILD_TYPE} STREQUAL "SERIAL" ) 
  set(ENV{FC} /shared/apps/rhel8/opt/rocm-afar-5891/bin/amdflang-new)
  set(CMAKE_Fortran_COMPILER ${FC} )
endif()     

set(FFTW_PATH "/shared/midgard/home/cedrick_eci/")
set(FFTW_LIB  "-L${FFTW_PATH}/lib -lfftw3")
# set(NCDF_LIB   "-lnetcdff") 
set(FFTW_INCLUDE_DIR "${FFTW_PATH}/include/")
set(LIBS             "${FFTW_LIB}")

#not used at the moment
#set(DRAGONEGG_FLAGS "-finline-aggressive -fslp-vectorize  -fmerge-all-constants") #  -mmadd4 -mfp64 -enable-strided-vectorization")

set(USER_Fortran_FLAGS         "-cpp ${USER_accel_FLAGS} ${USER_omp_FLAGS} -I${FFTW_INCLUDE_DIR}") #-fallow-argument-mismatch from gnu-version10
set(USER_Fortran_FLAGS_RELEASE "-O3") #these will be ignored:  -fprefetch-loop-arrays --param prefetch-latency=300") 
set(USER_Fortran_FLAGS_DEBUG   "-O0 -g -debug -ffpe-trap=all") 

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

add_definitions(-DUSE_FFTW ) # -DUSE_NETCDF) # -DUSE_BLAS -DUSE_MKL)
add_definitions(-DNO_ASSUMED_RANKS )


set(GNU_SED "gsed")

if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  add_definitions(-DUSE_PSFFT)

  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qopenmp")
  set(PSFFT_COMPILE_FLAGS "-qopenmp")

endif()
