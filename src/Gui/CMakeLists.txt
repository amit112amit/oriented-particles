project(Gui)

find_package(Qt5 COMPONENTS Core Gui Widgets)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SRCS LiveSimulationWindow.cxx LiveSimulation.cxx
    RunLiveSimulation.cxx)

qt5_wrap_cpp(HDR_MOC LiveSimulationWindow.h LiveSimulation.h)
qt5_wrap_ui(UI_MOC LiveSimulationWindow.ui)

add_executable(opsLive ${SRCS} ${HDR_MOC} ${UI_MOC})

# Link VTK, Qt and QWT libraries and include their headers
target_link_libraries(opsLive Qt5::Core Qt5::Gui
    Qt5::Widgets ${QWT_LIBRARIES} ${VTK_LIBRARIES})
target_include_directories(opsLive PUBLIC ${QWT_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Include OPS directories
set(INCLD Bodies Model Solvers Helpers)
foreach(DIR ${INCLD})
        target_link_libraries(opsLive ${DIR})
        target_include_directories(opsLive PRIVATE "${PROJECT_SOURCE_DIR}/../${DIR}")
endforeach(DIR)
