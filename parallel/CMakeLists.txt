FILE(GLOB sourcefiles "../parallel/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(parallel STATIC ${sourcefiles})
add_dependencies(parallel modules)
