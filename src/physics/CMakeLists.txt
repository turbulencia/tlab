FILE(GLOB sourcefiles "*.f90")
add_library(physics STATIC ${sourcefiles})
add_dependencies(physics mappings operators thermodynamics utils base)