project(src)

set(SUBD Bodies Model Solvers Drivers)

foreach(CURRD ${SUBD})
    add_subdirectory(${CURRD})
endforeach(CURRD)
