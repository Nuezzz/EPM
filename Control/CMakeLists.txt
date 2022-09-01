set(
    CONTROL_SRCS
    src/openmp.c
    src/timer.c
    src/error.c
    src/memory.c
    src/util.c
    src/vectormap.c
)

add_library(
    control
    ${CONTROL_SRCS}
)

target_link_libraries(control m)

target_include_directories(control PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")