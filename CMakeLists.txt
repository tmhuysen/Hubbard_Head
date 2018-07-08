cmake_minimum_required(VERSION 3.8)
project(Hubbard_Head VERSION 0.0.1
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# Use level-2 compiler optimization, since Armadillo really benefits from this
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Set CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/SetCMakeEnvironment.cmake)

# Find all required packages
include(${CMAKE_SOURCE_DIR}/cmake/FindPackages.cmake)

# Add the sources for the library
add_subdirectory(${PROJECT_SOURCE_FOLDER})

# Set up testing
enable_testing()
add_subdirectory(${PROJECT_TESTS_FOLDER})
