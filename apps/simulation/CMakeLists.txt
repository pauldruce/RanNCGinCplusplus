add_executable(Simulation Simulation.cpp SimulationData.hpp)

target_link_libraries(Simulation ${ARMADILLO_LIBRARIES})
target_link_libraries(Simulation RanNCG)

# We need this directory, and users of our library will need it too
target_include_directories(Simulation
                           PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/)
