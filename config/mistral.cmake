# Mistral cluster (DKRZ)

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" )
   set(ENV{FC} mpif90)
   set(CMAKE_Fortran_COMPILER mpif90) 
   set(USER_Fortran_FLAGS          " -fpp -nbs -convert little_endian -save-temps -xHost ") 
   set(USER_Fortran_FLAGS_RELEASE  " -O3")

   add_definitions(-DUSE_FFTW -DUSE_MPI -DUSE_MPI_IO )

   set(CMAKE_BUILD_TYPE RELEASE) 

else() # compiler for serial build
   set(ENV{FC} ifort)
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

set(FFTW_INCLUDE_DIR   "/sw/rhel6-x64/numerics/fftw-3.3.4-bullxmpi-gcc48/include")
set(FFTW_LIB           "/sw/rhel6-x64/numerics/fftw-3.3.4-bullxmpi-gcc48/lib/libfftw3.a")
set(LIBS ${FFTW_LIB} )

add_definitions(-DRESTRICTKEYWORD=__restrict__)
