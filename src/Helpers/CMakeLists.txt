project(Helpers)

# Add new Helper source file here
set(SRCS Kabsch.cxx Delaunay3DSurf.cxx)

add_library(Helpers ${SRCS})
target_include_directories(Helpers PUBLIC ${VTK_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS} )
