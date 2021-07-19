FILE(GLOB sourcefiles "../utils/*.f90")
add_library(utils STATIC ${sourcefiles})
add_dependencies(utils modules)
