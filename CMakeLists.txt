cmake_minimum_required(VERSION 3.15)

project(
  Edit-Distance
  VERSION 1.0
  DESCRIPTION "Edit Distance"
  HOMEPAGE_URL "https://github.com/ucrparlay/Edit-Distance"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_THREAD_LIBS_INIT "-pthread")

add_subdirectory(parlaylib)

option(DEBUG "Debug mode" OFF)

add_compile_options(-pthread)
if(NOT DEBUG)
  add_compile_options(-Wall -Wextra -O3 -mcx16 -march=native)
endif()

add_library(edit_distance_dp edit_distance_dp.cpp)
target_link_libraries(edit_distance_dp PUBLIC parlay)

add_library(edit_distance_parallel edit_distance_parallel.cpp)
target_link_libraries(edit_distance_parallel PUBLIC parlay)

add_executable(test_framework test_framework.cpp)
target_link_libraries(test_framework PUBLIC edit_distance_dp edit_distance_parallel)

add_library(suffix_array_sequential suffix_array_sequential.cpp)

add_executable(suffix_array_test suffix_array_test.cpp)
target_link_libraries(suffix_array_test PUBLIC suffix_array_sequential parlay)
