cmake_minimum_required(VERSION 3.16)
project(heat2d C)

set(CMAKE_C_STANDARD 99)

add_executable(lead_pdc_project cenatMD.c)
target_link_libraries(lead_pdc_project m)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(lead_pdc_project ${MPI_C_LIBRARIES})