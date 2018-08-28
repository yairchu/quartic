cmake_minimum_required(VERSION 3.0)
project(quartic_tests)

if(MSVC)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wno-missing-braces -std=c++11")
  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  endif()
endif()

add_subdirectory("../rapidcheck" build)

add_executable(quartic_tests_cpp
  quartic.c
  quartic_tests.cpp)
target_link_libraries(quartic_tests_cpp rapidcheck)
