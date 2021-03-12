FILE(GLOB sourcefiles "../filters/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(filters STATIC ${sourcefiles})
add_dependencies(filters modules)
