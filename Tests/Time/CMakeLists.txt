cmake_minimum_required(VERSION 2.6)
project (Tests)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

find_package(Armadillo REQUIRED)
find_package(yaml-cpp REQUIRED)

set(YAMLCPP_LIBRARY /usr/local/lib/libyaml-cpp.a)
set(INTEL_LIBRARY /opt/intel/mkl/lib/intel64/libmkl_rt.so)

add_executable(CholTime CholTime.cpp)
target_link_libraries(CholTime ${ARMADILLO_LIBRARIES} ${INTEL_LIBRARY} tbb)

