set(
    SRCS
    src/lapack_eigen.c
)

add_library(
    lapack_solver
    ${SRCS}
)

target_link_libraries(lapack_solver control ${MKL_CORE_LIB} ${MKL_GNU_LIB} ${MKL_ILP64_LIB})

target_include_directories(lapack_solver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" ${LAPACK_INCLUDE})
