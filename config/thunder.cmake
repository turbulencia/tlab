# Thunder cluster (MPIM) 

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build; assuming intel13 	  
   set(ENV{FC} mpif90)
   set(CMAKE_Fortran_COMPILER mpif90) 
   set(USER_Fortran_FLAGS          " -fpp -nbs -convert little_endian -save-temps -xHost ") 
   set(USER_Fortran_FLAGS_RELEASE  " -axavx,SSE4.2,SSE3,SSE2 -simd -vec-threshold0 -unroll-aggressive -opt-prefetch -O3")

   add_definitions(-DUSE_FFTW -DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL)

   set(CMAKE_BUILD_TYPE RELEASE) 

else() # compiler for serial build
   set(ENV{FC} /opt/local/bin/ifort)
   set(CMAKE_Fortran_COMPILER ifort) 
   set(USER_Fortran_FLAGS          "-fpp -DUSE_FFTW -nbs -save-temps -xHost")  

   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -convert big_endian")

     set(CMAKE_BUILD_TYPE RELEASE) 
     
   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" ) 
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -convert little_endian") 

     set(CMAKE_BUILD_TYPE RELEASE) 
     
   else()	
     set(USER_Fortran_FLAGS_DEBUG    "-Og -ggdb3 -Wall -Wno-unknown-pragmas")

     add_definitions(-D_DEBUG)
     
     set(CMAKE_BUILD_TYPE DEBUG) 

   endif() 

endif()     

set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/scratch/mpi/mh0738/cedrick/include")
set(FFTW_LIB           "/scratch/mpi/mh0738/cedrick/lib/libfftw3.a")
set(LIBS ${FFTW_LIB} )

add_definitions(-DRESTRICTKEYWORD=__restrict__)
