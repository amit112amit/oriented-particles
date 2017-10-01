cmake_minimum_required(VERSION 3.0.2)

project(OPS CXX Fortran)

cmake_policy(SET CMP0053 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(Eigen REQUIRED)
find_package(VTK REQUIRED)

include(${VTK_USE_FILE})

add_subdirectory(src)
