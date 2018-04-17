project(src)

set(SUBD Bodies Model Solvers Drivers Helpers Gui)

foreach(CURRD ${SUBD})
    add_subdirectory(${CURRD})
endforeach(CURRD)
