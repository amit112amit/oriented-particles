project(Helpers CXX Fortran)

# Add new Helper source file here
set(SRCS Kabsch.cxx Delaunay3DSurf.cxx ReadInputFile.cxx
    InterpolationFromData.cxx SimulationState.cxx StereoGraphicProjection.cxx)

add_library(Helpers STATIC ${SRCS})
target_link_libraries(Helpers PRIVATE Eigen3::Eigen vtkFiltersCore vtkFiltersExtraction vtkFiltersGeometry vtkIOLegacy)
