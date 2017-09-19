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

find_package(MKL)



add_executable(SELQR SELQR.cpp)
target_link_libraries(SELQR ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(QRSELQR QRSELQR.cpp)
target_link_libraries(QRSELQR ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(Comparation Comparation.cpp)
target_link_libraries(Comparation ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(QRSELQR_samplings QRSELQR_samplings.cpp)
target_link_libraries(QRSELQR_samplings ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${INTEL_LIBRARY} tbb)

add_executable(Quadroptor_simple Quadroptor_simple.cpp)

add_executable(Quadroptor Quadrotor.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS})
    target_link_libraries(Quadroptor_simple ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${MKL_LIBRARIES} pthread  tbb)
    target_link_libraries(Quadroptor ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} ${MKL_LIBRARIES} pthread  tbb)
else()
    target_link_libraries(Quadroptor_simple ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} tbb)
    target_link_libraries(Quadroptor ${ARMADILLO_LIBRARIES} ${YAMLCPP_LIBRARY} tbb)
endif()
