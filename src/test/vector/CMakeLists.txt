set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/.././bin/test/vector)

add_executable(vec_test vec_test.f90)
add_executable(v_test_2 vec_test_2.f90)

# Better Unit tests
add_executable(v_constructor_test v_constructor_test.f90)
add_executable(v_assignment_test v_assignment_test.f90)
add_executable(v_n_assigment_test v_n_assigment_test.f90)
add_executable(v_n_constructor_test v_n_constructor_test.f90)
add_executable(v_operator_test v_operator_test.f90)


target_link_libraries(vec_test PRIVATE formath)
target_link_libraries(v_test_2 PRIVATE formath)
target_link_libraries(v_constructor_test PRIVATE formath)
target_link_libraries(v_assignment_test PRIVATE formath)
target_link_libraries(v_n_assigment_test PRIVATE formath)
target_link_libraries(v_n_constructor_test PRIVATE formath)
target_link_libraries(v_operator_test PRIVATE formath)

add_test(NAME vector_test COMMAND vec_test)
add_test(NAME vector_fun COMMAND v_test_2)
add_test(NAME vector_constructor COMMAND v_constructor_test)
add_test(NAME vector_assignment COMMAND v_assignment_test)
add_test(NAME vector_n_assignment COMMAND v_n_assigment_test)
add_test(NAME vector_n_construction COMMAND v_n_constructor_test)
add_test(NAME vector_operator COMMAND v_operator_test)

add_custom_command(TARGET v_constructor_test COMMAND v_constructor_test COMMENT "Running constructor tests")
add_custom_command(TARGET v_assignment_test COMMAND v_assignment_test COMMENT "Running assigment tests")
add_custom_command(TARGET v_n_assigment_test COMMAND v_n_assigment_test COMMENT "Running n assigment tests")
add_custom_command(TARGET v_n_constructor_test COMMAND v_n_constructor_test COMMENT "Running n construction tests")
add_custom_command(TARGET v_operator_test COMMAND v_operator_test COMMENT "Running operator test")


