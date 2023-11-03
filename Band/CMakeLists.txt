set(
    SRCS
    src/epm_input.c
    src/epm.c
)

add_library(
    epm
    ${SRCS}
)

target_link_libraries(epm control reader lattice ham)

target_include_directories(epm PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")