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

option(DEBUG "" OFF)
option(CILKPLUS "" OFF)
option(OPENCILK "" OFF)
option(SERIAL "" OFF)
option(STDALLOC "" OFF)

if(CILKPLUS)
  add_compile_options(-DPARLAY_CILKPLUS -DCILK -fcilkplus)
elseif(OPENCILK)
  add_compile_options(-DPARLAY_OPENCILK -DCILK -fopencilk)
elseif(SERIAL)
  add_compile_options(-DPARLAY_SEQUENTIAL)
else()
  add_compile_options(-pthread)
endif()

if(DEBUG)
  set(CMAKE_BUILD_TYPE "Debug")
  add_compile_options(-Og -mcx16 -DDEBUG -g -Wall -Wextra)
else()
  set(CMAKE_BUILD_TYPE "Release")
  add_compile_options( -O3 -mcx16 -march=native)
endif()

if(STDALLOC)
  add_compile_options(-DPARLAY_USE_STD_ALLOC)
endif()

add_subdirectory(parlaylib)

add_library(edit_distance INTERFACE)
target_include_directories(edit_distance INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(edit_distance INTERFACE edit_distance_parallel suffix_array_sequential)

add_executable(test_framework test_framework.cpp)
target_link_libraries(test_framework PUBLIC edit_distance_parallel)

add_executable(test_framework_real test_framework_real.cpp)
target_link_libraries(test_framework_real PUBLIC edit_distance_parallel)

add_library(edit_distance_parallel edit_distance_parallel.cpp)
target_link_libraries(edit_distance_parallel PUBLIC parlay)

add_executable(suffix_array_test suffix_array_test.cpp)
target_link_libraries(suffix_array_test PUBLIC suffix_array_sequential parlay)

add_library(suffix_array_sequential suffix_array_sequential.cpp)
