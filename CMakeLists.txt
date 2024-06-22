cmake_minimum_required(VERSION 3.10)
project(X17_n_TOF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

set(ROOT_DIR ${PROJECT_SOURCE_DIR})

include_directories(
    ${ROOT_DIR}/code/sipm
    ${ROOT_DIR}/code/urwell
    ${ROOT_INCLUDE_DIRS}
)

file(GLOB SIPM_SOURCES ${ROOT_DIR}/code/sipm/*.cpp)
file(GLOB URWELL_SOURCES ${ROOT_DIR}/code/urwell/*.cpp)
file(GLOB CUBE_SOURCES ${ROOT_DIR}/code/cube/*.cpp)
set(UTILS_SOURCE ${ROOT_DIR}/code/utils.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -O3")

add_library(sipm_library ${SIPM_SOURCES})
add_library(urwell_library ${URWELL_SOURCES})
add_library(cube_library ${CUBE_SOURCES})
add_library(utils_library ${UTILS_SOURCE})

add_executable(sipm_analysis ${ROOT_DIR}/code/sipm/sipm_analysis.cpp)
target_link_libraries(sipm_analysis sipm_library cube_library utils_library ${ROOT_LIBRARIES})

add_executable(urwell_analysis ${ROOT_DIR}/code/urwell/urwell_analysis.cpp)
target_link_libraries(urwell_analysis urwell_library cube_library utils_library ${ROOT_LIBRARIES})
