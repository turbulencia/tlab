FILE(GLOB sourcefiles "*.f90")
add_library(statistics STATIC ${sourcefiles})
add_dependencies(statistics physics operators thermodynamics utils ibm base)