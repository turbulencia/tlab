FILE(GLOB sourcefiles "../mappings/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(mappings STATIC ${sourcefiles})
add_dependencies(mappings modules)
