FILE(GLOB sourcefiles "../io/*.f90")
add_library(io STATIC ${sourcefiles})
add_dependencies(io modules)
