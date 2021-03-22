FILE(GLOB sourcefiles "../mpi/*.f90")
add_library(mpi STATIC ${sourcefiles})
add_dependencies(mpi modules)
