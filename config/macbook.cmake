# MacBook

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build	  
   set(ENV{FC} /opt/local/bin/mpif90-gcc48)   
   set(CMAKE_Fortran_COMPILER /opt/local/bin/mpif90-gcc48)
   set(USER_Fortran_FLAGS          "-ffree-form -ffree-line-length-2048 -fno-automatic ")
   set(USER_Fortran_FLAGS_RELEASE  "-O3 -cpp -arch x86_64 -ffast-math -ffinite-math-only -funroll-loops -mtune=native ")   

   add_definitions(-DUSE_PSFFT -DUSE_MPI -DUSE_ALLTOALL -DUSE_MPI_IO)

   set(CMAKE_BUILD_TYPE RELEASE) 

else() # compiler for serial build
   set(ENV{FC} /opt/local/bin/gfortran-mp-4.8)
   set(USER_Fortran_FLAGS          "-cpp -ffree-form -ffree-line-length-2048 -fno-automatic")  

   add_definitions(-DUSE_PSFFT -DUSE_FFTW)

   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -fconvert=big-endian    -mtune=native -ffast-math -ffinite-math-only -funroll-loops")

     set(CMAKE_BUILD_TYPE RELEASE) 
     
   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" ) 
     set(USER_Fortran_FLAGS_RELEASE  "-O3 -fconvert=little-endian -mtune=native -ffast-math -ffinite-math-only -funroll-loops") 

     set(CMAKE_BUILD_TYPE RELEASE) 
     
   else()	
     set(USER_Fortran_FLAGS_DEBUG    "-Og -ggdb3 -Wall -Wno-unknown-pragmas")

     add_definitions(-D_DEBUG)

     set(CMAKE_BUILD_TYPE DEBUG) 

   endif() 

endif()     

set(GNU_SED "gsed")

set(FFTW_INCLUDE_DIR   "/opt/local/include")
set(FFTW_LIB           "/opt/local/lib/libfftw3.a")
#set(FFTW_INCLUDE_DIR   "/usr/local/fftw_intel14/include")
#set(FFTW_LIB           "/usr/local/fftw_intel14/lib/libfftw3.a")
set(LIBS ${FFTW_LIB} )

add_definitions(-DRESTRICTKEYWORD=__restrict__)
