project(Bodies)

# Append name of the source file to SRCS
set(SRCS OPSBody.cxx BrownianBody.cxx ViscosityBody.cxx ALConstraint.cxx OPSMesh.cxx)

add_library(Bodies ${SRCS})
target_include_directories(Bodies PUBLIC ${VTK_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
target_include_directories(Bodies PUBLIC "${PROJECT_SOURCE_DIR}/../Helpers")
