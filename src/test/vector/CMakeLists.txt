set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/.././bin/test/vector)

set(VECTOR_TEST_FILES

    TestVectorAssignment
    TestVectorAssignmentSpeed
    TestVectorConstructor
    TestVectorConstructorSpeed
    TestVectorOperator
    TestVectorOperatorSpeed
    TestVectorSubroutineOperator
    TestVectorSubroutineOperatorSpeed
    TestVectorMisc

)

foreach(tfile ${VECTOR_TEST_FILES})

    add_executable(${tfile} ${tfile}.f90)
    add_test(NAME ${tfile} COMMAND ${tfile})
    target_link_libraries(${tfile} formath)

endforeach(tfile)



