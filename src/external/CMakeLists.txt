FILE(GLOB sourcefiles "fitpack/*.f90")
add_library(fitpack STATIC ${sourcefiles})
add_dependencies(fitpack modules)

FILE(GLOB sourcefiles "quadpack/*.f90")
add_library(quadpack STATIC ${sourcefiles})
add_dependencies(quadpack modules)
