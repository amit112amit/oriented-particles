project(Model)

add_library(Model Model.cxx)
target_include_directories(Model PUBLIC ${EIGEN_INCLUDE_DIRS})
target_include_directories(Model PUBLIC "${PROJECT_SOURCE_DIR}/../Bodies/")
