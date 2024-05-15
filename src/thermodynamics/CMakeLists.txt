FILE(GLOB sourcefiles "*.f90")
add_library(thermodynamics STATIC ${sourcefiles})
add_dependencies(thermodynamics io modules)
