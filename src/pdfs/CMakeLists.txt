FILE(GLOB sourcefiles "../pdfs/*.f90")
add_library(pdfs STATIC ${sourcefiles})
add_dependencies(pdfs modules)
