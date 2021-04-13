FILE(GLOB sourcefiles "../operators/*.f90")
add_library(operators STATIC ${sourcefiles})
add_dependencies(operators modules)
