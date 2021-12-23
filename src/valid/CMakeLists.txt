set(LEVEL1_LIBS utils modules)

FILE(GLOB sourcefiles "partial/*.f90")
add_executable(vpartial.x ${sourcefiles})
target_link_libraries(vpartial.x operators io fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "integral/*.f90")
add_executable(vintegral.x ${sourcefiles})
target_link_libraries(vintegral.x operators io fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "poisson/*.f90")
add_executable(vpoisson.x ${sourcefiles})
target_link_libraries(vpoisson.x mappings operators io thermodynamics filters fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "burgers/*.f90")
add_executable(vburgers.x ${sourcefiles})
target_link_libraries(vburgers.x mappings operators io thermodynamics fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "interpol/*.f90")
add_executable(vinterpol.x ${sourcefiles})
target_link_libraries(vinterpol.x operators io fdm ${LEVEL1_LIBS} ${LIBS})