project(PostProcess CXX Fortran)

add_library(shtoolswrapper STATIC SHTools_Wrapper.f95)
target_include_directories( shtoolswrapper PRIVATE ${SHTOOLS_INCLUDE_DIR} )
target_link_libraries( shtoolswrapper PRIVATE ${SHTOOLS_LIBRARY}
    ${FFTW3_LIBRARY})

add_executable(postProc PostProcess.cxx)
target_link_libraries(postProc Eigen3::Eigen ${CGAL_LIBRARIES} shtoolswrapper
    ${Boost_LIBRARIES} ${VTK_LIBRARIES})
target_include_directories(postProc PRIVATE ${CGAL_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS})

add_dependencies( postProc shtoolswrapper )

add_executable(plotRow PlotRow.cxx)
target_link_libraries(plotRow ${VTK_LIBRARIES} Eigen3::Eigen Helpers)
target_include_directories(plotRow PRIVATE ${VTK_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/../Helpers)
