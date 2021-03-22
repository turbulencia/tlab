FILE(GLOB sourcefiles "../thermodynamics/*.f90")
add_library(thermodynamics STATIC ${sourcefiles})
add_dependencies(thermodynamics modules)
