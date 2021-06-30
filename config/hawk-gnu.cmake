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
  set(ENV{FC} mpif90 )
  set(CMAKE_Fortran_COMPILER mpif90)
  #set(CMAKE_EXE_LINKER_FLAGS "-flto") 
  #set(CMAKE_AR  "gcc-ar")
  #set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
  #set(CMAKE_RANLIB "")

  
  #set(ENV{FC} gfortran)
  #set(CMAKE_Fortran_COMPILER gfortran) 
endif()     

set(USER_Fortran_FLAGS         "-cpp -ffree-form -ffree-line-length-2048 -fno-automatic") #-fallow-argument-mismatch from gnu-version10
set(USER_Fortran_FLAGS_RELEASE "-march=znver2 -mtune=znver2 -O3 -ffinite-math-only -fprefetch-loop-arrays --param prefetch-latency=300")
# Do not use -funroll-all-loops / -funroll-loops 
set(USER_Fortran_FLAGS_DEBUG   "-g -traceback -debug all -ffpe-trap=all") 

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

add_definitions(-DUSE_FFTW -DHLRS_HAWK -DUSE_NETCDF) # -DUSE_BLAS -DUSE_MKL) 
set(FFTW_LIB   "-lfftw3")
set(MKL_LIB    "")
set(NCDF_LIB   "-lnetcdff") 

set(LIBS             "${NCDF_LIB} ${FFTW_LIB} ${MKL_LIB} -l:libamdlibm.a -lm")

set(GNU_SED "gsed")

if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  add_definitions(-DUSE_PSFFT)

  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qopenmp")
  set(PSFFT_COMPILE_FLAGS "-qopenmp")

endif()
