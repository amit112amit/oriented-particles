project(src)

set(SUBD Bodies Model Solvers Drivers Helpers Constraints)

foreach(CURRD ${SUBD})
    add_subdirectory(${CURRD})
endforeach(CURRD)
