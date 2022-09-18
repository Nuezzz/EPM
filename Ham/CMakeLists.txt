set(
    SRCS
    src/ham_local.c
    src/ham_spin.c
)

add_library(
    ham
    ${SRCS}
)

target_link_libraries(ham lattice control lapack_solver)

target_include_directories(ham PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")