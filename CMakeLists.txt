cmake_minimum_required(VERSION 3.0.2)

project(OPS CXX Fortran)

cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0054 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets)

find_package(VTK
 COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersExtraction
    FiltersGeneral
    FiltersGeometry
    FiltersSources
    IOLegacy
    RenderingCore
    RenderingOpenGL2
    RenderingFreeType
    InteractionStyle
    RenderingGL2PSOpenGL2
)

find_package(CGAL REQUIRED)
include( ${CGAL_USE_FILE} )

find_package(Boost REQUIRED COMPONENTS "filesystem")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_library( FFTW3_LIBRARY "fftw3" HINTS ${FFTW3_ROOT} )

find_library( SHTOOLS_LIBRARY "SHTOOLS" HINTS ${SHTOOLS_ROOT} )
if(${SHTOOLS_LIBRARY} STREQUAL "SHTOOLS-NOTFOUND")
        message(" SHTOOLS_LIBRARY not found in default locations. Set SHTOOLS_ROOT to specify custom location." )
endif( ${SHTOOLS_LIBRARY} STREQUAL "SHTOOLS-NOTFOUND" )

find_path( SHTOOLS_INCLUDE_DIR "shtools.mod" PATHS ${SHTOOLS_ROOT} PATH_SUFFIXES "include")
if(${SHTOOLS_INCLUDE_DIR} STREQUAL "SHTOOLS_INCLUDE_DIR-NOTFOUND")
        message(" File shtools.mod was not found in default locations. Set SHTOOLS_ROOT to specify custom location." )
endif( ${SHTOOLS_INCLUDE_DIR} STREQUAL "SHTOOLS_INCLUDE_DIR-NOTFOUND" )

#find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# Set path for qwt
if(CMAKE_HOST_UNIX)
	# qwt has different locations on different linux distributions
	find_path(QWT_INCLUDE_DIRS "qwt_series_data.h"
		  HINTS "/usr/include/" "/usr/local/include"
		  PATH_SUFFIXES "qwt" "qt5/qwt6"
		  DOC "qwt header files location.")
	find_library(QWT_LIBRARIES NAMES "qwt" "qwt-qt5"
		  HINTS "/usr/lib/" "/usr/local/lib"
		  DOC "qwt shared library.")
endif(CMAKE_HOST_UNIX)
if(CMAKE_HOST_WIN32)
	set(QWT_INCLUDE_DIRS "C:/msys64/mingw64/include/qwt" CACHE PATH "Set Path for qwt header files.")
	set(QWT_LIBRARIES "C:/msys64/mingw64/lib/libqwt.dll.a" CACHE STRING "Name of qwt library.")
endif(CMAKE_HOST_WIN32)

# We will add two resource files for opsLive gui demo
install(FILES resources/T7.vtk resources/T7_OPS_Asphericity.dat TYPE BIN)

add_subdirectory(src)
