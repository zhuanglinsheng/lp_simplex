add_executable(test_simple_examples test_simple_examples.c)
target_compile_definitions(test_simple_examples PRIVATE
	SIMPLE_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/data/simple_examples"
)
target_link_libraries(test_simple_examples PRIVATE lp_simplex)
add_test(NAME test_simple_examples COMMAND test_simple_examples)

add_executable(test_degeneracy test_degeneracy.c)
target_include_directories(test_degeneracy PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_degeneracy PRIVATE lp_simplex)
add_test(NAME test_degeneracy COMMAND test_degeneracy)

add_executable(test_presolve test_presolve.c)
target_include_directories(test_presolve PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_presolve PRIVATE lp_simplex)
add_test(NAME test_presolve COMMAND test_presolve)
