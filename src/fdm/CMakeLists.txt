FILE(GLOB sourcefiles "../fdm/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(fdm STATIC ${sourcefiles})
add_dependencies(fdm modules)
