# Optionally glob, but only for CMake 3.12 or later:
file(GLOB HEADER_LIST CONFIGURE_DEPENDS
     "${RandomNCG_SOURCE_DIR}/include/RanNCG/*.h")
# set(HEADER_LIST "${ModernCMakeExample_SOURCE_DIR}/include/modern/lib.hpp")

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
target_include_directories(RanNCG
                           PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/)

# set_target_properties(RanNCG PROPERTIES CXX_STANDARD 17)
target_include_directories(RanNCG PUBLIC ${RandomNCG_SOURCE_DIR}/dependencies)
