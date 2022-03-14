add_executable(
  DefaultActionSimulations
  DefaultActionSimulations.cpp
  ${RandomNCG_SOURCE_DIR}/include/RanNCG/SimulationData.hpp)

target_link_libraries(DefaultActionSimulations RanNCG)

# We need this directory, and users of our library will need it too
target_include_directories(
  DefaultActionSimulations PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/
                                  ${HDF5_INCLUDE_DIR})
