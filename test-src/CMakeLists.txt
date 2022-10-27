# Example test.
add_executable(example example.cpp)
target_link_libraries(example
        ${MPFR_LIBRARIES})
