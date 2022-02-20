add_executable(Simulation Simulation.cpp ../../include/RanNCG/SimulationData.hpp)

#find_package(HDF5 REQUIRED COMPONENTS CXX)

target_link_libraries(Simulation RanNCG)

# We need this directory, and users of our library will need it too
target_include_directories(
  Simulation PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/ ${HDF5_INCLUDE_DIR})
