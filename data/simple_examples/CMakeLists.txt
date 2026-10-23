set(simple_example_names
	production_planning
	mixed_inequalities
	equality_system
	free_variables
	free_and_lower_bounds
	variable_bounds
	ranged_rows
	degenerate_pivot
	energy_system
)

foreach(name IN LISTS simple_example_names)
	add_executable(example_${name}_manual ${name}_1.c)
	add_executable(example_${name}_mps ${name}_2.c)
	foreach(mode IN ITEMS manual mps)
		target_link_libraries(example_${name}_${mode} PRIVATE lp_simplex)
		if(NOT MSVC)
			target_link_libraries(example_${name}_${mode} PRIVATE m)
		endif()
	endforeach()
	if(BUILD_TESTING)
		add_test(NAME example_${name}_manual
			COMMAND example_${name}_manual)
		add_test(NAME example_${name}_mps
			COMMAND example_${name}_mps ${CMAKE_CURRENT_SOURCE_DIR}/${name}.mps)
	endif()
endforeach()
