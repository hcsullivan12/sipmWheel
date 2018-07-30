cmake_minimum_required(VERSION 3.5)
project(sipmwheel)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT COMPONENTS REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

set(INCLUDE_DIRECTORIES ${ROOT_INCLUDE_DIR})
include_directories( ${INCLUDE_DIRECTORIES})
set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR})
link_directories( ${LINK_DIRECTORIES})

link_libraries(stdc++fs)

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.h)
add_executable(sipmwheel main.cxx ${sources} ${headers})

target_link_libraries(sipmwheel Spectrum)
target_link_libraries(sipmwheel ${ROOT_LIBRARIES})
