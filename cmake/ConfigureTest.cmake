# In this CMake file, we will include header files and link to libraries for a given test source

# for each test:
# ... add the boost headers ...
target_include_directories(${TEST_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# ... add this project's library ...
target_include_directories(${TEST_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
target_link_libraries(${TEST_NAME} PRIVATE ${LIBRARY_NAME})

# ... add Armadillo ...
target_include_directories(${TEST_NAME} PUBLIC ${ARMADILLO_INCUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC ${ARMADILLO_LIBRARIES})

# ... add Eigen3 ...
target_include_directories(${TEST_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC ${EIGEN3_LIBRARIES})
