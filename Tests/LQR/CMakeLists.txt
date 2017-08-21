cmake_minimum_required(VERSION 2.6)
project (Tests)

SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
SET (CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib /opt/intel/mkl/lib/intel64/)

find_package(Armadillo REQUIRED)
find_package(yaml-cpp)

link_directories(/usr/local/lib)

set(YAMLCPP_LIBRARY yaml-cpp)
set(INTEL_LIBRARY /opt/intel/mkl/lib/intel64/libmkl_rt.so)


add_executable(SELQR SELQR.cpp)
target_link_libraries(SELQR ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(QRSELQR QRSELQR.cpp)
target_link_libraries(QRSELQR ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(Comparation Comparation.cpp)
target_link_libraries(Comparation ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(QRSELQR_samplings QRSELQR_samplings.cpp)
target_link_libraries(QRSELQR_samplings ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(Quadroptor_simple Quadroptor_simple.cpp)
target_link_libraries(Quadroptor_simple ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)
