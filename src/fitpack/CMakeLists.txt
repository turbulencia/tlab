FILE(GLOB sourcefiles "../fitpack/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(fitpack STATIC ${sourcefiles})
add_dependencies(fitpack modules)
