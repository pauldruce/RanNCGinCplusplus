set(HEADER_LIST
    ${RandomNCG_SOURCE_DIR}/include/RanNCG/Clifford.h
    ${RandomNCG_SOURCE_DIR}/include/RanNCG/DiracMatrices.h
    ${RandomNCG_SOURCE_DIR}/include/RanNCG/DiracOperator.h
    ${RandomNCG_SOURCE_DIR}/include/RanNCG/Pauli.h
    ${RandomNCG_SOURCE_DIR}/include/RanNCG/RandomMatrixUtils.h
    ${RandomNCG_SOURCE_DIR}/include/RanNCG/Simulation.h
    ${RandomNCG_SOURCE_DIR}/include/RanNCG/SimulationData.hpp)

add_library(RanNCG)

target_sources(
  RanNCG
  PRIVATE Clifford.cpp
          DiracMatrices.cpp
          DiracOperator.cpp
          Pauli.cpp
          Simulation.cpp
          RandomMatrixUtils.cpp
          README.md
  PUBLIC ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(RanNCG PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/
                                         ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(RanNCG PRIVATE ${ARMADILLO_LIBRARIES})

target_include_directories(
  RanNCG PUBLIC ${RandomNCG_SOURCE_DIR}/dependencies/cppitertools)

#target_include_directories(RanNCG
#                           PUBLIC ${RandomNCG_SOURCE_DIR}/dependencies/doctest)
