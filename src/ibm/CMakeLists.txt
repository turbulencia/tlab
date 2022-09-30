FILE(GLOB sourcefiles "../ibm/*.f90")
add_library(ibm STATIC ${sourcefiles})
add_dependencies(ibm io utils modules)