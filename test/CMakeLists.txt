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

add_executable(test_pan test_pan.c)
target_link_libraries(test_pan PRIVATE lp_simplex)
add_test(NAME test_pan COMMAND test_pan)

add_executable(test_pan_dynamic_sparse_qr test_pan_dynamic_sparse_qr.c)
target_include_directories(test_pan_dynamic_sparse_qr PRIVATE
	${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_dynamic_sparse_qr PRIVATE lp_simplex)
add_test(NAME test_pan_dynamic_sparse_qr COMMAND test_pan_dynamic_sparse_qr)

if(LP_SIMPLEX_USE_SPQR AND SPQR_INCLUDE_DIR AND SPQR_LIBRARY)
	add_executable(test_pan_multifrontal_qr test_pan_multifrontal_qr.c)
	target_include_directories(test_pan_multifrontal_qr PRIVATE
		${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
	target_include_directories(test_pan_multifrontal_qr SYSTEM PRIVATE
		${SPQR_INCLUDE_DIR})
	target_link_libraries(test_pan_multifrontal_qr PRIVATE lp_simplex
		${CHOLMOD_LIBRARY} ${SUITESPARSE_CONFIG_LIBRARY})
	add_test(NAME test_pan_multifrontal_qr COMMAND test_pan_multifrontal_qr)
endif()

add_executable(test_presolve test_presolve.c)
target_include_directories(test_presolve PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_presolve PRIVATE lp_simplex)
add_test(NAME test_presolve COMMAND test_presolve)
