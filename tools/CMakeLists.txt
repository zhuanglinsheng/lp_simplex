add_executable(test_netlib test_netlib.c)
target_compile_definitions(test_netlib PRIVATE
	NETLIB_REFERENCE_CSV="${PROJECT_SOURCE_DIR}/data/netlib/feasible_gurobi_1e-8.csv"
)
target_link_libraries(test_netlib PRIVATE lp_simplex)

if(BUILD_TESTING AND EXISTS "${PROJECT_SOURCE_DIR}/data/netlib/feasible/afiro.mps")
	add_test(
		NAME netlib_afiro_smoke
		COMMAND test_netlib
			${PROJECT_SOURCE_DIR}/data/netlib/feasible/afiro.mps
	)
	add_test(
		NAME netlib_afiro_dual_revised
		COMMAND test_netlib --algorithm dual-revised
			${PROJECT_SOURCE_DIR}/data/netlib/feasible/afiro.mps
	)
	add_test(
		NAME netlib_afiro_pan_bda
		COMMAND test_netlib --algorithm pan-bda
			${PROJECT_SOURCE_DIR}/data/netlib/feasible/afiro.mps
	)
endif()

if(BUILD_TESTING AND EXISTS "${PROJECT_SOURCE_DIR}/data/netlib/feasible/degen2.mps")
	add_test(
		NAME netlib_degen2_dual_revised
		COMMAND test_netlib --algorithm dual-revised
			${PROJECT_SOURCE_DIR}/data/netlib/feasible/degen2.mps
	)
endif()

# Numerical regressions which exercise the raw-model recovery and alternate
# dual-crash paths.  Keep these as individual tests so a future failure names
# the affected model directly.
foreach(model IN ITEMS 25fv47 agg2 bnl1 etamacro grow22 maros scfxm1 share1b)
	if(BUILD_TESTING AND
	   EXISTS "${PROJECT_SOURCE_DIR}/data/netlib/feasible/${model}.mps")
		add_test(
			NAME netlib_${model}_dual_revised
			COMMAND test_netlib --algorithm dual-revised
				${PROJECT_SOURCE_DIR}/data/netlib/feasible/${model}.mps
		)
	endif()
endforeach()
