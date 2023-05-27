FILE(GLOB sourcefiles "*.f90")
add_library(operators STATIC ${sourcefiles})
add_dependencies(operators ibm filters fdm modules)