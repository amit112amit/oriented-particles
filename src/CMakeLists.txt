project(src)

set(SUBD Bodies Model Solvers Drivers Helpers Gui PostProcess)

foreach(CURRD ${SUBD})
    add_subdirectory(${CURRD})
endforeach(CURRD)
