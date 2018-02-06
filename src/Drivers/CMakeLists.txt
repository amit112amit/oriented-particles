project(Drivers)

#Set the include directory names which are under ../../src/
set(INCLD Bodies Model Solvers Helpers)

#Add a new add_executable for each driver
add_executable(opsConsistency TestOPSConsistency.cxx)
add_executable(heatOPS BrownDynOPS.cxx)
add_executable(asphericity Asphericity.cxx)
add_executable(transition CrumpleOrMelt.cxx)
#add_executable(plate PeriodicPlateMelting.cxx)

# Add HDF5 libraries and headers
#find_package(HDF5 REQUIRED COMPONENTS CXX)

#Add name of your executable target to the DRIVERS variable
#set(DRIVERS opsConsistency asphericity heatOPS transition plate)
set(DRIVERS opsConsistency asphericity heatOPS transition)

#For each driver link libraries and set the include paths
foreach(DRIVER ${DRIVERS})    
    target_link_libraries(${DRIVER} ${VTK_LIBRARIES} ${HDF5_LIBRARIES})
    target_include_directories(${DRIVER} PRIVATE ${VTK_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
    foreach(DIR ${INCLD})
        target_link_libraries(${DRIVER} ${DIR})
        target_include_directories(${DRIVER} PRIVATE "${PROJECT_SOURCE_DIR}/../${DIR}")
    endforeach(DIR)
endforeach(DRIVER)
