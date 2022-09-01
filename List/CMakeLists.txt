set(
    SRCS
    src/list.c
)

add_library(
    list
    ${SRCS}
)

target_link_libraries(list control)

target_include_directories(list PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")