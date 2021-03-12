FILE(GLOB sourcefiles "../quadpack/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(quadpack STATIC ${sourcefiles})
add_dependencies(quadpack modules)
