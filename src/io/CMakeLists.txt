FILE(GLOB sourcefiles "../io/*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
MESSAGE( STATUS "Compile defs contain: " ${defs} )
if ( "USE_NETCDF" IN_LIST defs )
else()
FILE(GLOB to_remove "../io/io_averages.f90") # get the complete path
list(REMOVE_ITEM sourcefiles ${to_remove})
endif()

add_library(io STATIC ${sourcefiles})
add_dependencies(io modules)
