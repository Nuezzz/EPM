set(
    SRCS
    src/reader.c
    src/reader_utils.c
)

add_library(
    reader
    ${SRCS}
)

target_link_libraries(reader list)

target_include_directories(reader PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")