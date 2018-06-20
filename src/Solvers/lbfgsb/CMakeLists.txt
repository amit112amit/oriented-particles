project(lbfgsb Fortran)

IF(WIN32 AND MSVC)
	IF(NOT ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel"))
		MESSAGE("Please set MINGW_GFORTRAN to mingw-gfortran.exe location.")
		set (MINGW_GFORTRAN "C:/Software/mingw-w64/mingw64/bin/gfortran.exe" CACHE FILEPATH "Please provide full path to mingw-gfortran.exe")
	ENDIF()
ENDIF(WIN32 AND MSVC)

add_library(lbfgsb STATIC lbfgsb-routines.f)
