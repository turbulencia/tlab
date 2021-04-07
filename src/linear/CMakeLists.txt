FILE(GLOB sourcefiles "../linear/*.f90")
add_library(linear STATIC ${sourcefiles})
add_dependencies(linear modules)
