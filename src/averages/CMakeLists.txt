FILE(GLOB sourcefiles "../averages/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(averages STATIC ${sourcefiles})
add_dependencies(averages modules)
