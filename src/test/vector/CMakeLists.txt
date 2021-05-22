set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/.././bin/test/vector)

add_executable(vec_test vec_test.f90)
add_executable(v_test_2 vec_test_2.f90)

target_link_libraries(vec_test PRIVATE formath)
target_link_libraries(v_test_2 PRIVATE formath)

add_test(NAME vector_assignment COMMAND vec_test)
add_test(NAME vector_fun COMMAND v_test_2)


