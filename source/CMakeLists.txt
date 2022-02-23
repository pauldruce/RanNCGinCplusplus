add_library(
  RanNCG
  Clifford.cpp
  DiracMatrices.cpp
  DiracOperator.cpp
  Pauli.cpp
  Simulation.cpp
  RandomMatrixUtils.cpp
  README.md
  ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(RanNCG PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/
                                         ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(RanNCG PRIVATE ${ARMADILLO_LIBRARIES})

target_include_directories(RanNCG PUBLIC ${RandomNCG_SOURCE_DIR}/dependencies)
