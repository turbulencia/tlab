FILE(GLOB sourcefiles "../thermodynamics/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(thermodynamics STATIC ${sourcefiles})
add_dependencies(thermodynamics modules)
