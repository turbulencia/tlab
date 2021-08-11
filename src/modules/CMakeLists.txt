FILE(GLOB sourcefiles "*.f90")

# retrieve the git hash from the current commit
find_package(Git)
if(GIT_EXECUTABLE)
  exec_program(git ${CMAKE_CURRENT_SOURCE_DIR} ARGS "describe --dirty" OUTPUT_VARIABLE GITHASH RETURN_VALUE GITRETURN)
  if(NOT GITRETURN MATCHES 0)
    set(GITHASH "not available")
  endif()
else()
  set(GITHASH "not available")
endif()

message(STATUS "Git hash " ${GITHASH})

# send a precompiler statement replacing the git hash
add_definitions(-DGITHASH="${GITHASH}")

get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_MPI" IN_LIST defs )
else()
  FILE(GLOB to_remove "tlab_mpi_vars.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
  FILE(GLOB to_remove "tlab_mpi_procs.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()

add_library(modules ${sourcefiles})
