add_executable(RanNCGTest  CliffordTests.cpp
        DiracOperatorTests.cpp
        Type11SimulationTests.cpp)

target_include_directories(
  RanNCGTest
  PUBLIC ${RandomNCG_SOURCE_DIR}/include
         ${RandomNCG_SOURCE_DIR}/apps/simulation
         ${DOCTEST_INCLUDE_DIR})
target_link_libraries(RanNCGTest RanNCG)
