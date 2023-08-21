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
  # include(CTest)
endif()

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# set(DEPENDENCIES_DIR "dependencies") set(DEPENDENCIES
# ${DEPENDENCIES_DIR}/cppitertools/combinations.hpp)
# include_directories(${DEPENDENCIES_DIR})

# Get DOCTEST
include(ExternalProject)
find_package(Git REQUIRED)

#ExternalProject_Add(
#  doctest
#  PREFIX ${CMAKE_BINARY_DIR}/doctest
#  GIT_REPOSITORY https://github.com/doctest/doctest.git
#  TIMEOUT 10
#  UPDATE_COMMAND ${GIT_EXECUTABLE} pull
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND ""
#  INSTALL_COMMAND ""
#  LOG_DOWNLOAD ON)
#
## Expose required variable (DOCTEST_INCLUDE_DIR) to parent scope
#ExternalProject_Get_Property(doctest source_dir)
#set( DOCTEST_INCLUDE_DIR
#    ${source_dir}/doctest
#    CACHE INTERNAL "Path to include folder for doctest")

### Get DOCTEST
include(FetchContent)
FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        GIT_TAG v2.4.8)
FetchContent_MakeAvailable(doctest)
set(DOCTEST_INCLUDE_DIR
        "${doctest_SOURCE_DIR}/"
        CACHE INTERNAL "Path to include folder for doctest")

include_directories(DOCTEST_INCLUDE_DIR)

# Build project library
add_subdirectory(source)

# Build executables that use library
add_subdirectory(apps/simulation)
add_subdirectory(apps/data_analysis)

# Add tests if building tests.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(tests)
endif()

# install(TARGETS Type11Simulation DESTINATION bin)
