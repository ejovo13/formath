set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/.././bin/test/matrix)

set(MATRIX_TEST_FILES 
    TestMatrixAssignment 
    TestMatrixConstructors 
    TestMatrixOperators 
    TestMatrixSubroutineOperators
    TestMatrixQRFactorization    
)

foreach(tfile ${MATRIX_TEST_FILES})

    add_executable(${tfile} ${tfile}.f90)
    add_test(NAME ${tfile} COMMAND ${tfile})
    target_link_libraries(${tfile} Formath)

endforeach(tfile)
