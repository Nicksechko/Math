cmake_minimum_required(VERSION 3.14)
project(Math)

set(CMAKE_CXX_STANDARD 17)

add_executable(Math main.cpp fraction.cpp fraction.h matrix.h matrix.cpp linear_system.cpp lup_matrix.cpp options.cpp extended_matrix.cpp permutation.cpp)