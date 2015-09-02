if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE BIG)
endif()
 
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build; not developed
   message( FATAL_ERROR "Only BIG mode developed. CMake will exit." )

else() # compiler for serial build
   set(ENV{FC} gfortran)
   set(CMAKE_Fortran_COMPILER gfortran) 
   set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic")

   add_definitions(-DUSE_FFTW)

   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
     set(USER_Fortran_FLAGS_RELEASE "-fconvert=big-endian -O3 -ffast-math -mtune=native -march=native")
     
     set(CMAKE_BUILD_TYPE RELEASE) 
     
   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" ) 
     set(USER_Fortran_FLAGS_RELEASE "-fconvert=little-endian -O3 -ffast-math -mtune=native -march=native")
     
     set(CMAKE_BUILD_TYPE RELEASE) 
     
   else()
     set(USER_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -Wno-unknown-pragmas -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal")
     
     add_definitions(-D_DEBUG)
     
     set(CMAKE_BUILD_TYPE DEBUG) 
     
   endif() 

endif()     

set(FFTW_INCLUDE_DIR   "/usr/local/include")
set(FFTW_LIB           "/usr/local/lib/libfftw3.a")
set(LIBS ${FFTW_LIB})
