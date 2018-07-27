project(Model)

add_library(Model Model.cxx)
target_link_libraries(Model Eigen3::Eigen)
target_include_directories(Model PUBLIC "${PROJECT_SOURCE_DIR}/../Bodies/")
