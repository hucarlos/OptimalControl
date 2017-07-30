cmake_minimum_required(VERSION 2.6)
project (Tests)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

find_package(Armadillo REQUIRED)
find_package(yaml-cpp REQUIRED)

set(YAMLCPP_LIBRARY /usr/local/lib/libyaml-cpp.a)

add_executable(Quadratic Quadratic.cpp)
target_link_libraries(Quadratic ${ARMADILLO_LIBRARIES})

add_executable(Obstacles Obstacles.cpp)
target_link_libraries(Obstacles ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY})

add_executable(SystemCost SystemCost)
target_link_libraries(SystemCost ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY})
