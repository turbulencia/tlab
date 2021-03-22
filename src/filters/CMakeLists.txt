FILE(GLOB sourcefiles "../filters/*.f90")
add_library(filters STATIC ${sourcefiles})
add_dependencies(filters modules)
