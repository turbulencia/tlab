FILE(GLOB sourcefiles "../mappings/*.f90")
add_library(mappings STATIC ${sourcefiles})
add_dependencies(mappings modules)
