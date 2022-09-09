set(
    SRCS
    src/lattice_input.c
    src/lattice_neighbor.c
)

add_library(
    lattice
    ${SRCS}
)

target_link_libraries(lattice reader control)

target_include_directories(lattice PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")