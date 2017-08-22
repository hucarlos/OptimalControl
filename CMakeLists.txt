project(LQRAlgorithms)
cmake_minimum_required(VERSION 2.8)


SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
 
find_package(Armadillo REQUIRED)
find_package(yaml-cpp REQUIRED)


set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(LOCAL_INCLUDE_DIR /usr/local/include)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)

  
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
  
include_directories("${PROJECT_INCLUDE_DIR}" ${ARMADILLO_INCLUDE_DIRS} ${yaml-cpp_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${LOCAL_INCLUDE_DIR}")

link_directories(/usr/local/lib)


add_subdirectory(Tests)






