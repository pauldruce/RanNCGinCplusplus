cmake_minimum_required(VERSION 3.7...3.22)
if(${CMAKE_VERSION} VERSION_LESS 3.22)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.22)
endif()

project(
  RandomNCG
  VERSION 1.0
  DESCRIPTION
    "A library and application to run a monte carlo simulation of finite non-commutative geometries."
  LANGUAGES CXX)

# Only do these if this is the main project, and not if it is included through
# add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)

  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

  # Testing only available if this is the main app Note this needs to be done in
  # the main CMakeLists since it calls enable_testing, which must be in the main
  # CMakeLists. include(CTest)
endif()

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# set(DEPENDENCIES_DIR "dependencies") set(DEPENDENCIES
# ${DEPENDENCIES_DIR}/cppitertools/combinations.hpp)
# include_directories(${DEPENDENCIES_DIR})

add_subdirectory(source)
add_subdirectory(apps/simulation)
add_subdirectory(apps/data_analysis)

install(TARGETS Type11Simulation DESTINATION bin)
