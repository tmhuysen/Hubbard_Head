# Configuration file for the "Hubbard_Head" package
# It defines the following variables:
#   hubbard_head_INCLUDE_DIRS   - header include directories
#   hubbard_head_LIBRARIES      - library to link against


# Specify the include directory
set(hubbard_head_INCLUDE_DIRS /usr/local/hubbard_head/include)

# Import the exported targets
include(/usr/local/hubbard_head/cmake/Hubbard_HeadTargets.cmake)

# Specify the library value
set(hubbard_head_LIBRARIES Hubbard_Head)

# Output the library version
message(STATUS "hubbard_head version: 0.0.1")
