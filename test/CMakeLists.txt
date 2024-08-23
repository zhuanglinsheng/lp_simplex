include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(test_lp_simplex_1 test_lp_simplex_1.c)
target_link_libraries(test_lp_simplex_1 lp_simplex)
add_test(
	NAME test_lp_simplex_1
	COMMAND test_lp_simplex_1
)

add_executable(test_lp_simplex_2 test_lp_simplex_2.c)
target_link_libraries(test_lp_simplex_2 lp_simplex)
add_test(
	NAME test_lp_simplex_2
	COMMAND test_lp_simplex_2
)

add_executable(test_lp_simplex_3 test_lp_simplex_3.c)
target_link_libraries(test_lp_simplex_3 lp_simplex)
add_test(
	NAME test_lp_simplex_3
	COMMAND test_lp_simplex_3
)

add_executable(test_lp_simplex_4 test_lp_simplex_4.c)
target_link_libraries(test_lp_simplex_4 lp_simplex)
add_test(
	NAME test_lp_simplex_4
	COMMAND test_lp_simplex_4
)

add_executable(test_lp_simplex_5 test_lp_simplex_5.c)
target_link_libraries(test_lp_simplex_5 lp_simplex)
add_test(
	NAME test_lp_simplex_5
	COMMAND test_lp_simplex_5
)

add_executable(test_lp_simplex_6 test_lp_simplex_6.c)
target_link_libraries(test_lp_simplex_6 lp_simplex)
add_test(
	NAME test_lp_simplex_6
	COMMAND test_lp_simplex_6
)

add_executable(test_lp_simplex_7 test_lp_simplex_7.c)
target_link_libraries(test_lp_simplex_7 lp_simplex)
add_test(
	NAME test_lp_simplex_7
	COMMAND test_lp_simplex_7
)
