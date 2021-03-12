FILE(GLOB sourcefiles "../linear/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(linear STATIC ${sourcefiles})
add_dependencies(linear modules)
