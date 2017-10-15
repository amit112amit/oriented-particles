project(Drivers)

#Set the include directory names which are under ../../src/
set(INCLD Bodies Model Solvers Helpers)

#Add a new add_executable for each driver
add_executable(opsConsistency TestOPSConsistency.cxx)
add_executable(heatOPS BrownDynOPS.cxx)
add_executable(freeDiff FreeDiffusion.cxx)
add_executable(constrained BrownDynOPSConstraints.cxx)

#Add name of your executable target to the DRIVERS variable
set(DRIVERS opsConsistency heatOPS freeDiff constrained)

#For each driver link libraries and set the include paths
foreach(DRIVER ${DRIVERS})    
    target_link_libraries(${DRIVER} ${VTK_LIBRARIES})
    target_include_directories(${DRIVER} PRIVATE ${VTK_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
    foreach(DIR ${INCLD})
        target_link_libraries(${DRIVER} ${DIR})
        target_include_directories(${DRIVER} PRIVATE "${PROJECT_SOURCE_DIR}/../${DIR}")
    endforeach(DIR)
endforeach(DRIVER)

