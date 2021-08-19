FILE(GLOB sourcefiles "*.f90")
add_library(fdm STATIC ${sourcefiles})
add_dependencies(fdm modules)
