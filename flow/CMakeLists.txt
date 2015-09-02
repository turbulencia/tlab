FILE(GLOB sourcefiles "../flow/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(flow STATIC ${sourcefiles})
add_dependencies(flow modules)
