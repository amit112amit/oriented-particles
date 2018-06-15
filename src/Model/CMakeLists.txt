project(Model)

add_library(Model Model.cxx OPSModel.cxx)
target_include_directories(Model PUBLIC ${EIGEN_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS} ${CGAL_LIBRARIES})
target_include_directories(Model PUBLIC "${PROJECT_SOURCE_DIR}/../Bodies/" "${PROJECT_SOURCE_DIR}/../Helpers/")
