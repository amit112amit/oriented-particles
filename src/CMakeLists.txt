project(src)

set(SUBD Bodies Model Solvers Drivers Helpers PostProcess)

foreach(CURRD ${SUBD})
    add_subdirectory(${CURRD})
endforeach(CURRD)
