set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/.././bin/test)

add_executable(vec_test vec_test.f90)
add_executable(gs_test gram_schmidt_test.f90)
add_executable(gs_test_2 gs_test_2.f90)
add_executable(v_test_2 vec_test_2.f90)

target_compile_options(gs_test PRIVATE -Ofast)

target_link_libraries(vec_test PRIVATE formath)
target_link_libraries(gs_test PRIVATE formath)
target_link_libraries(gs_test_2 PRIVATE formath)
target_link_libraries(v_test_2 PRIVATE formath)