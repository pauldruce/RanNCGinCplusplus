add_executable(RanNCGTests CliffordTests.cpp)
# Includes
target_include_directories(
  RanNCGTests PUBLIC ${DOCTEST_INCLUDE_DIR}
                     ${RandomNCG_SOURCE_DIR}/include/RanNCG)

# Libary linking
target_link_libraries(RanNCGTests RanNCG)
