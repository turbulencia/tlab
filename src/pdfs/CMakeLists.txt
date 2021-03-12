FILE(GLOB sourcefiles "../pdfs/*.f90")
include_directories(${INCLUDE_DIRS} "../include")

add_library(pdfs STATIC ${sourcefiles})
add_dependencies(pdfs modules)
