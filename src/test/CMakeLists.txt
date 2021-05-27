add_executable(gs_test gram_schmidt_test.f90)
add_executable(gs_test_2 gs_test_2.f90)

target_compile_options(gs_test PRIVATE -Ofast)

add_subdirectory(vector)
add_subdirectory(matrix)

target_link_libraries(gs_test PRIVATE Formath)
target_link_libraries(gs_test_2 PRIVATE Formath)