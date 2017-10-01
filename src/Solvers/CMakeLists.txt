project(Solvers)

include(CMakeAddFortranSubdirectory)
cmake_add_fortran_subdirectory( lbfgsb 		
		PROJECT lbfgsb
		ARCHIVE_DIR ../
		RUNTIME_DIR ../
		LIBRARIES lbfgsb
		LINK_LIBRARIES
		NO_EXTERNAL_INSTALL
		)

set(SOLVERS LBFGSBWrapper.cxx)

add_library(Solvers ${SOLVERS})
target_include_directories(Solvers PUBLIC "${EIGEN_INCLUDE_DIRS}")
target_include_directories(Solvers PUBLIC "${PROJECT_SOURCE_DIR}/../Model")
target_include_directories(Solvers PUBLIC "${PROJECT_SOURCE_DIR}/../Bodies")

