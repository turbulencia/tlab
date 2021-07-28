FILE(GLOB sourcefiles "*.f90")
add_library(mappings STATIC ${sourcefiles})
add_dependencies(mappings thermodynamics utils modules)
