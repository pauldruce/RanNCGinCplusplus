add_executable(AverageEigenvalues average_eigenvalues.cpp)

target_link_libraries(AverageEigenvalues PRIVATE RanNCG)

# We need this directory, and users of our library will need it too
target_include_directories(AverageEigenvalues
                           PUBLIC ${RandomNCG_SOURCE_DIR}/include/RanNCG/)
