cmake_minimum_required(VERSION 3.16)
project(parsum C)

set(CMAKE_C_STANDARD 99)

add_executable(parsum parsum.c)
target_link_libraries(parsum m)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(parsum ${MPI_C_LIBRARIES})