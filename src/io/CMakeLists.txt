FILE(GLOB sourcefiles "*.f90")
add_library(io STATIC ${sourcefiles})
add_dependencies(io utils base)
