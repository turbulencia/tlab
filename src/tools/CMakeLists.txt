set(LEVEL1_LIBS utils base)
if (PSFFT_LIB)
set(LEVEL1_LIBS ${LEVEL1_LIBS} ${PSFFT_LIB} ${FFTWTHREADS_LIB})
endif()

# Simulation
FILE(GLOB sourcefiles "dns/*.f90")

if (PSFFT_LIB)
else()
FILE(GLOB to_remove "dns/rhs_global_incompressible_nbc.f90") # get the complete path
list(REMOVE_ITEM sourcefiles ${to_remove})
endif()

add_executable(dns.x ${sourcefiles})
target_link_libraries(dns.x statistics physics mappings operators particles filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_source_files_properties(dns/rhs_global_incompressible_nbc.f90 PROPERTIES COMPILE_FLAGS "${PSFFT_COMPILE_FLAGS}")
set_target_properties(dns.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

# Initialization
FILE(GLOB sourcefiles "initialize/grid/*.f90")
add_executable(inigrid.x ${sourcefiles})
target_link_libraries(inigrid.x operators fdm io ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "initialize/rand/*.f90")
add_executable(inirand.x ${sourcefiles})
target_link_libraries(inirand.x statistics physics mappings operators filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(inirand.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

FILE(GLOB sourcefiles "initialize/scal/*.f90")
add_executable(iniscal.x ${sourcefiles})
target_link_libraries(iniscal.x statistics physics mappings operators filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(iniscal.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

FILE(GLOB sourcefiles "initialize/flow/*.f90")
add_executable(iniflow.x ${sourcefiles})
target_link_libraries(iniflow.x statistics physics mappings operators filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(iniflow.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

# Postprocessing
FILE(GLOB sourcefiles "plot/visuals*.f90")
add_executable(visuals.x ${sourcefiles})
target_link_libraries(visuals.x statistics physics mappings operators particles filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(visuals.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

FILE(GLOB sourcefiles "statistics/averages*.f90")
add_executable(averages.x ${sourcefiles})
target_link_libraries(averages.x statistics physics mappings operators particles filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(averages.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

FILE(GLOB sourcefiles "statistics/spectra*.f90")
add_executable(spectra.x ${sourcefiles})
target_link_libraries(spectra.x statistics physics mappings operators filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(spectra.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

FILE(GLOB sourcefiles "statistics/pdfs*.f90")
add_executable(pdfs.x ${sourcefiles})
target_link_libraries(pdfs.x statistics physics mappings operators filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(pdfs.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

# Transformation
FILE(GLOB sourcefiles "transform/transgrid.f90")
add_executable(transgrid.x ${sourcefiles})
target_link_libraries(transgrid.x operators io ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "transform/transfields.f90")
add_executable(transfields.x ${sourcefiles})
target_link_libraries(transfields.x statistics physics mappings operators filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(transfields.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

# Particles
FILE(GLOB sourcefiles "initialize/part/*.f90")
add_executable(inipart.x ${sourcefiles})
target_link_libraries(inipart.x statistics physics mappings operators particles filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(inipart.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

# FILE(GLOB sourcefiles "initialize/part/particle_ini_trajec.f90")
# add_executable(l_ini_trajec.x ${sourcefiles})
# target_link_libraries(l_ini_trajec.x statistics physics mappings operators particles fdm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_ini_trajec.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "particles/particle_build_pdf.f90")
# add_executable(l_pdf.x ${sourcefiles})
# target_link_libraries(l_pdf.x statistics physics mappings operators particles filters fdm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_pdf.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "particles/particle_trajec.f90")
# add_executable(l_trajec.x ${sourcefiles})
# target_link_libraries(l_trajec.x statistics physics mappings operators particles filters fdm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_trajec.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "particles/particle_pos_trajec*.f90")
# add_executable(l_pos_trajec.x ${sourcefiles})
# target_link_libraries(l_pos_trajec.x statistics physics mappings operators particles filters fdm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_pos_trajec.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# Structure analysis
FILE(GLOB sourcefiles "structure/apriori.f90")
add_executable(apriori.x ${sourcefiles})
target_link_libraries(apriori.x statistics physics mappings operators filters fdm ibm thermodynamics io ${LEVEL1_LIBS} ${LIBS})
if (PSFFT_LIB)
set_target_properties(apriori.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
endif()

# Cloud
FILE(GLOB sourcefiles "cloud/state.f90")
add_executable(state.x ${sourcefiles})
target_link_libraries(state.x thermodynamics io ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "cloud/smooth.f90")
add_executable(smooth.x ${sourcefiles})
target_link_libraries(smooth.x thermodynamics io ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "cloud/saturation.f90")
add_executable(saturation.x ${sourcefiles})
target_link_libraries(saturation.x thermodynamics io ${LEVEL1_LIBS} ${LIBS})
