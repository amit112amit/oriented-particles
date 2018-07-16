cmake_minimum_required(VERSION 3.0.2)

project(OPS CXX Fortran)

cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0054 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen REQUIRED)
find_package(VTK REQUIRED)
find_package(Boost REQUIRED)
find_package(CGAL REQUIRED)
#find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()
include(${VTK_USE_FILE})
include(${CGAL_USE_FILE})

add_subdirectory(src)
