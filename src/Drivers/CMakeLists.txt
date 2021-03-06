project(Drivers)

#Set the include directory names which are under ../../src/
set(INCLD Bodies Model Solvers Helpers)

#Add a new add_executable for each driver
add_executable(opsConsistency TestOPSConsistency.cxx)
add_executable(heatOPS BrownDynOPS.cxx)
add_executable(opsModel HeatOPSModel.cxx)
add_executable(asphericity Asphericity.cxx)
add_executable(inertia Inertia.cxx)
add_executable(plate PlateMelting.cxx)
add_executable(printDoF PrintDOF.cxx)

# Add HDF5 libraries and headers
#find_package(HDF5 REQUIRED COMPONENTS CXX)

# Add name of your executable target to the DRIVERS variable
set(DRIVERS opsConsistency asphericity heatOPS inertia plate opsModel printDoF)

# For each driver link libraries and set the include paths
foreach(DRIVER ${DRIVERS})
	target_link_libraries(${DRIVER} PRIVATE ${HDF5_LIBRARIES} Eigen3::Eigen
		VTK::IOLegacy)
	target_include_directories(${DRIVER} PRIVATE ${HDF5_INCLUDE_DIRS})
	foreach(DIR ${INCLD})
		target_link_libraries(${DRIVER} PRIVATE ${DIR})
		target_include_directories(${DRIVER} PRIVATE "${PROJECT_SOURCE_DIR}/../${DIR}")
	endforeach(DIR)
endforeach(DRIVER)

# Install the executables on `make install` or equivalent
install(TARGETS ${DRIVERS} RUNTIME OPTIONAL)
