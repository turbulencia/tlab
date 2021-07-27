FILE(GLOB sourcefiles "*.f90")
add_library(utils STATIC ${sourcefiles})
add_dependencies(utils modules)
