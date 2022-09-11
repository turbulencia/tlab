FILE(GLOB sourcefiles "*.f90")

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
else()
FILE(GLOB to_remove "particle_sort.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
FILE(GLOB to_remove "particle_send_recv.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
endif()

add_library(particles STATIC ${sourcefiles})
add_dependencies(particles thermodynamics io modules)
