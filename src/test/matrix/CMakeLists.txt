set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/.././bin/test/matrix)

add_executable(m_asgn_test mat_asgn_test.f90)
add_executable(m_ctor_test mat_ctor_test.f90)
add_executable(m_oper_test mat_oper_test.f90)
add_executable(m_subop_test mat_subop_test.f90)

target_link_libraries(m_asgn_test formath)
target_link_libraries(m_ctor_test formath)
target_link_libraries(m_oper_test formath)
target_link_libraries(m_subop_test formath)

add_test(NAME "Matrix_Assignment_Test" COMMAND m_asgn_test)
add_test(NAME "Matrix_Constructor_Test" COMMAND m_ctor_test)
add_test(NAME "Matrix_Operator_Test" COMMAND m_oper_test)
add_test(NAME "Matrix_Sub_Operator_Test" COMMAND m_subop_test)