FILE(GLOB sourcefiles "../io/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(io STATIC ${sourcefiles})
add_dependencies(io modules)
