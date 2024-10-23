FILE(GLOB sourcefiles "*.f90")
add_library(filters STATIC ${sourcefiles})
add_dependencies(filters base)
