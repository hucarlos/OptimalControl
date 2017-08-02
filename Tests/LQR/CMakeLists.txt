cmake_minimum_required(VERSION 2.6)
project (Tests)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

find_package(Armadillo REQUIRED)
find_package(yaml-cpp REQUIRED)

FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()

link_directories(/usr/local/lib)

set(YAMLCPP_LIBRARY yaml-cpp)


add_executable(SELQR SELQR.cpp)
target_link_libraries(SELQR ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY})

add_executable(QRSELQR QRSELQR.cpp)
target_link_libraries(QRSELQR ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY})

add_executable(Comparation Comparation.cpp)
target_link_libraries(Comparation ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY})
