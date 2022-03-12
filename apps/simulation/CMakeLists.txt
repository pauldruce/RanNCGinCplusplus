add_executable(
  Type11Simulation Type11Simulation.cpp
                   ${RandomNCG_SOURCE_DIR}/include/RanNCG/SimulationData.hpp)

target_link_libraries(Type11Simulation RanNCG)

# We need this directory, and users of our library will need it too
target_include_directories(
  Type11Simulation PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/
                          ${HDF5_INCLUDE_DIR})
