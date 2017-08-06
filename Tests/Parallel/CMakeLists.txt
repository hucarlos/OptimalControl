cmake_minimum_required(VERSION 2.6)
project (Tests)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

set (PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_INCLUDE_DIR}")


add_executable(TBB1 TBB1.cpp)
target_link_libraries(TBB1 ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} tbb)

add_executable(Regression Regression.cpp)
target_link_libraries(Regression ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} tbb)
