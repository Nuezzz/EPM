add_library(
    dfise
    src/dfise.c
)

add_library(
    gmsh
    src/gmsh_convert.c
)

target_link_libraries(dfise geometry)
target_link_libraries(gmsh geometry)
