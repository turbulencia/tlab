set(LEVEL1_LIBS utils modules)

FILE(GLOB sourcefiles "partial/*.f90")
add_executable(vpartial.x ${sourcefiles})
target_link_libraries(vpartial.x operators ibm fitpack io thermodynamics fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "integral/*.f90")
add_executable(vintegral.x ${sourcefiles})
target_link_libraries(vintegral.x operators ibm fitpack io thermodynamics fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "poisson/*.f90")
add_executable(vpoisson.x ${sourcefiles})
target_link_libraries(vpoisson.x mappings operators ibm fitpack io thermodynamics filters fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "burgers/*.f90")
add_executable(vburgers.x ${sourcefiles})
target_link_libraries(vburgers.x mappings operators ibm fitpack io thermodynamics fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "spline/*.f90")
add_executable(vspline.x ${sourcefiles})
target_link_libraries(vspline.x operators ibm fitpack io fdm ${LEVEL1_LIBS} ${LIBS})
