########################################################
# Linux Specific Intel Compiler
# 
# -C         array index checking at run-time and more
# -fpe0      stop on floating-point exceptions
# -traceback give location of error at the end
# -p         compiles for profiling with gprof
#            To be set in compilation and linking
# -g         debug information
# -fast      Macro to increase performance, but 
#            numerical results are changed. See SUN block.
#            To be set in compilation and linking
# -warn interfaces to check arguments of routines
# -warn unused     to check declared variables not used
# -DUSE_FFTW Use FFTW library. 
#            Add path FFTWROOT as needed
# -openmp    Use OPENMP (compiler and linker flag)
#            Add -USE_OPENMP and remove -save
# -DUSE_MPI Use USE_MPI sections. Add also -DUSE_MPI_IO.
#            At MPI-M, use wrapper mpif90 instead of ifort
#            The mpi library at MPI-M works only in little-endian!
########################################################
if ( NOT BUILD_TYPE )
   message( WARNING "Setting BUILD_TYPE to default value." )
   set(BUILD_TYPE BIG)
endif()
 
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build; assuming intel13
   set(ENV{FC} mpif90)
   set(CMAKE_Fortran_COMPILER mpif90) 
   set(USER_Fortran_FLAGS          " -fpp -save -nbs -convert little_endian -xHost ") 
   set(USER_Fortran_FLAGS_RELEASE  " -O3 ")

   add_definitions(-DUSE_FFTW -DUSE_MPI -DUSE_MPI_IO)

   set(CMAKE_BUILD_TYPE RELEASE)

else() # compiler for serial build
   set(ENV{FC} ifort)
   set(CMAKE_Fortran_COMPILER ifort) 
   set(USER_Fortran_FLAGS             " -fpp -save -nbs -xHost -warn unused")  

   add_definitions(-DUSE_FFTW)

   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -convert big_endian")
     
     set(CMAKE_BUILD_TYPE RELEASE) 
     
   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" ) 
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -convert little_endian") 
     
     set(CMAKE_BUILD_TYPE RELEASE) 
     
   else()	
     #      set(USER_Fortran_FLAGS_DEBUG    "-Og -ggdb3 -Wall -Wno-unknown-pragmas")
     set(USER_Fortran_FLAGS_DEBUG    "-O0 -pg -debug full -traceback -fpe0 -save-temps -Wall")
     
     add_definitions(-D_DEBUG)
     
     set(CMAKE_BUILD_TYPE DEBUG) 
     
   endif() 

endif()     

set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/include")
set(FFTW_LIB           "/sw/squeeze-x64/numerics/fftw-3.3-openmp-gccsys/lib/libfftw3.a")
set(LIBS ${FFTW_LIB} )

add_definitions(-DRESTRICTKEYWORD=__restrict__)
