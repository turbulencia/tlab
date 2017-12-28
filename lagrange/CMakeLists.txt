FILE(GLOB sourcefiles "../lagrange/*.f90")

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
else()
FILE(GLOB to_remove "../lagrange/particle_sort.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
FILE(GLOB to_remove "../lagrange/particle_send_recv.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
FILE(GLOB to_remove "../lagrange/particle_to_field_send_recv_east.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
FILE(GLOB to_remove "../lagrange/particle_to_field_send_recv_north.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
endif()

include_directories(${INCLUDE_DIRS} "../include")

add_library(lagrange STATIC ${sourcefiles})
add_dependencies(lagrange modules)
