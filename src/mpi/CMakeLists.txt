FILE(GLOB sourcefiles "../mpi/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(mpi STATIC ${sourcefiles})
add_dependencies(mpi modules)
