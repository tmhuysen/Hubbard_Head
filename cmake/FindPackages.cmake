# In this CMake file, we will find all required packages


# Find the boost package - needed for unittests
find_package(Boost REQUIRED)

# Find Armadillo for linear algebra operations
find_package(Armadillo REQUIRED)

# Find Eigen3
find_package(Eigen3 REQUIRED)
