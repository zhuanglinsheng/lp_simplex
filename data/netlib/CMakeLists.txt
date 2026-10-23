# Dataset-specific references stay here; lp_solve_mps accepts explicit expectations.
if(NOT BUILD_TESTING)
    return()
endif()

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/feasible_gurobi_1e-8.csv" references)

function(add_netlib_test name model)
    set(mps "${CMAKE_CURRENT_SOURCE_DIR}/feasible/${model}.mps")
    if(NOT EXISTS "${mps}")
        return()
    endif()
    foreach(reference IN LISTS references)
        string(REPLACE "," ";" fields "${reference}")
        list(GET fields 0 reference_model)
        list(GET fields 1 reference_status)
        if(reference_model STREQUAL model AND reference_status STREQUAL "OPTIMAL")
            list(GET fields 3 objective)
            add_test(NAME ${name}
                COMMAND lp_solve_mps ${ARGN} --expected ${objective} "${mps}")
            return()
        endif()
    endforeach()
    message(FATAL_ERROR "Missing Netlib reference objective for ${model}")
endfunction()

add_netlib_test(netlib_afiro_smoke afiro)
add_netlib_test(netlib_afiro_dual_revised afiro --algorithm dual-revised)
add_netlib_test(netlib_afiro_pan_bda afiro --algorithm pan-bda)

# Keep numerical regressions individually named for diagnostics.
foreach(model IN ITEMS degen2 25fv47 agg2 bnl1 etamacro grow22 maros perold pilot.ja pilot.we pilot87 qap8 scfxm1 share1b tuff)
    add_netlib_test(netlib_${model}_dual_revised ${model} --algorithm dual-revised)
endforeach()
