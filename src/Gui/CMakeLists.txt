project(Gui)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SRCS1 LiveSimulationWindow.cxx LiveSimulation.cxx
	RunLiveSimulation.cxx)

set(SRCS2 RoughnessWidget.cxx RunRoughnessGUI.cxx)

qt5_wrap_cpp(HDR_MOC1 LiveSimulationWindow.h LiveSimulation.h)
qt5_wrap_ui(UI_MOC1 LiveSimulationWindow.ui)

qt5_wrap_cpp(HDR_MOC2 RoughnessWidget.h)
qt5_wrap_ui(UI_MOC2 RoughnessGUI.ui)

add_executable(opsLive ${SRCS1} ${HDR_MOC1} ${UI_MOC1})
add_executable(roughness ${SRCS2} ${HDR_MOC2} ${UI_MOC2})

# Link VTK, Qt and QWT libraries and include their headers
set(EXECS opsLive roughness)
# Include OPS directories
set(INCLD Bodies Model Solvers Helpers)
foreach(EXE ${EXECS})
	target_link_libraries(${EXE} PRIVATE ${QWT_LIBRARIES} Eigen3::Eigen
		VTK::FiltersSources VTK::GUISupportQt VTK::IOLegacy VTK::RenderingCore
		VTK::RenderingOpenGL2 VTK::RenderingFreeType VTK::InteractionStyle
		VTK::RenderingGL2PSOpenGL2)
	target_include_directories(${EXE} PRIVATE ${QWT_INCLUDE_DIRS}
		${Boost_INCLUDE_DIRS})
	foreach(DIR ${INCLD})
		target_link_libraries(${EXE} PRIVATE ${DIR})
		target_include_directories(${EXE} PRIVATE "${PROJECT_SOURCE_DIR}/../${DIR}")
	endforeach(DIR)
endforeach(EXE)

install(TARGETS opsLive roughness RUNTIME OPTIONAL)
