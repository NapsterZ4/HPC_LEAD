cmake_minimum_required(VERSION 3.16)
project(heat2d C)

set(CMAKE_C_STANDARD 99)

add_executable(heat2d heat2d.c)
target_link_libraries(heat2d m)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(heat2d ${MPI_C_LIBRARIES})