project(Constraints)

# Append name of the source file to SRCS
set(SRCS AugmentedLagrangian.cxx)

add_library(Constraints ${SRCS})
target_include_directories(Constraints PUBLIC ${EIGEN_INCLUDE_DIRS})
target_include_directories(Constraints PUBLIC "${PROJECT_SOURCE_DIR}/../Bodies")
