FILE(GLOB sourcefiles "../operators/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(operators STATIC ${sourcefiles})
add_dependencies(operators modules)
