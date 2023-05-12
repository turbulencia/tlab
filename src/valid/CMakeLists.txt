set(LEVEL1_LIBS utils modules)

FILE(GLOB sourcefiles "efilter/*.f90")
add_executable(vefilter.x efilter/vefilter.f90)
target_link_libraries(vefilter.x operators filters ibm fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "partial/*.f90")
add_executable(vpartial.x ${sourcefiles})
target_link_libraries(vpartial.x operators ibm io fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "integral/*.f90")
add_executable(vintegral.x ${sourcefiles})
target_link_libraries(vintegral.x operators ibm io fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "poisson/*.f90")
add_executable(vpoisson.x ${sourcefiles})
target_link_libraries(vpoisson.x mappings operators ibm io thermodynamics filters fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "helmholtz_fxz/*.f90")
add_executable(vhelmholtz_fxz.x ${sourcefiles})
target_link_libraries(vhelmholtz_fxz.x mappings operators ibm io thermodynamics filters fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "burgers/*.f90")
add_executable(vburgers.x ${sourcefiles})
target_link_libraries(vburgers.x mappings operators filters ibm io thermodynamics fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "interpol/*.f90")
add_executable(vinterpol.x ${sourcefiles})
target_link_libraries(vinterpol.x operators ibm io fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "interpartial/*.f90")
add_executable(vinterpartial.x ${sourcefiles})
target_link_libraries(vinterpartial.x mappings operators ibm io thermodynamics fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "cspline/*.f90")
add_executable(vcspline.x ${sourcefiles})
target_link_libraries(vcspline.x operators fdm ${LEVEL1_LIBS} ${LIBS})