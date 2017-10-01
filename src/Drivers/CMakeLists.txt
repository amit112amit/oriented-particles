project(Drivers)

#Set the include directory names which are under ../../src/
set(INCLD Bodies Model Solvers)

#Add a new add_executable for each driver
add_executable(thermal BrownDynDriver.cxx)

#Add name of your executable to the DRIVERS variable
set(DRIVERS thermal)

#For each driver link libraries and set the include paths
foreach(DRIVER ${DRIVERS})    
    target_link_libraries(${DRIVER} lbfgsb ${VTK_LIBRARIES})
    target_include_directories(${DRIVER} PRIVATE ${VTK_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
    foreach(DIR ${INCLD})
        target_link_libraries(${DRIVER} ${DIR})
        target_include_directories(${DRIVER} PRIVATE "${PROJECT_SOURCE_DIR}/../${DIR}")
    endforeach(DIR)
endforeach(DRIVER)

