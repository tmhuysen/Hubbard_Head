# In this CMake file, we will include the headers and link to the necessary libraries


# Include this project's headers
target_include_directories(${LIBRARY_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})

# Include the boost headers
target_include_directories(${LIBRARY_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# Include Armadillo
target_include_directories(${LIBRARY_NAME} PUBLIC ${ARMADILLO_INCUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC ${ARMADILLO_LIBRARIES})

# Include Eigen
target_include_directories(${LIBRARY_NAME} PUBLIC ${Eigen3_INCUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC ${Eigen3_LIBRARIES})
