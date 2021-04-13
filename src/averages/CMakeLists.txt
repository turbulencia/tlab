FILE(GLOB sourcefiles "../averages/*.f90")
add_library(averages STATIC ${sourcefiles})
add_dependencies(averages modules)
