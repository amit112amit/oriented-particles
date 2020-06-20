project(Bodies)

# Append name of the source file to SRCS
set(SRCS OPSBody.cxx BrownianBody.cxx ViscosityBody.cxx ALConstraint.cxx
        OPSMesh.cxx Morse2D.cxx Pressure.cxx OPSModel.cxx)

add_library(Bodies ${SRCS})
target_link_libraries(Bodies PRIVATE ${CGAL_LIBRARIES}
	VTK::FiltersExtraction VTK::FiltersGeometry VTK::IOLegacy Eigen3::Eigen)
	
target_include_directories(Bodies PUBLIC ${CGAL_INCLUDE_DIRS})
target_include_directories(Bodies PUBLIC "${PROJECT_SOURCE_DIR}/../Helpers")
