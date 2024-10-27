FILE(GLOB sourcefiles "*.f90")
add_library(physics STATIC ${sourcefiles})
add_dependencies(physics operators thermodynamics utils base)