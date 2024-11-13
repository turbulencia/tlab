FILE(GLOB sourcefiles "*.f90")
add_library(statistics STATIC ${sourcefiles})
add_dependencies(statistics physics mappings operators thermodynamics utils ibm base)