FILE(GLOB sourcefiles "../modules/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(modules ${sourcefiles})
