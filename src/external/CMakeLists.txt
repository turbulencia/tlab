FILE(GLOB sourcefiles "fitpack/*.f90")
add_library(fitpack STATIC ${sourcefiles})
add_dependencies(fitpack modules)

