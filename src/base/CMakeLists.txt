FILE(GLOB sourcefiles "*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_MPI" IN_LIST defs )
else()
  FILE(GLOB to_remove "tlab_mpi_vars.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
  FILE(GLOB to_remove "tlab_mpi_procs.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()
add_library(base ${sourcefiles})
