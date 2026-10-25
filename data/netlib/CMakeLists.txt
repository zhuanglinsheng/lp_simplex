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
        if(reference_model STREQUAL model AND reference_status STREQUAL
                                              "OPTIMAL")
            list(GET fields 3 objective)
            if("--experimental-environment" IN_LIST ARGN)
                add_test(
                    NAME ${name}
                    COMMAND
                        ${CMAKE_COMMAND} "-DDRIVER=$<TARGET_FILE:lp_solve_mps>"
                        "-DDRIVER_ARGUMENTS=${ARGN}" "-DEXPECTED=${objective}"
                        "-DMODEL=${mps}" -P
                        "${PROJECT_SOURCE_DIR}/test/verify_cli_configuration.cmake"
                )
            else()
                add_test(NAME ${name} COMMAND lp_solve_mps ${ARGN} --expected
                                              ${objective} "${mps}")
            endif()
            return()
        endif()
    endforeach()
    message(FATAL_ERROR "Missing Netlib reference objective for ${model}")
endfunction()

add_netlib_test(netlib_afiro_smoke afiro)
add_netlib_test(netlib_afiro_dual_revised afiro --algorithm dual-revised)
add_netlib_test(netlib_afiro_pan_bda afiro --algorithm pan-bda)
add_netlib_test(netlib_afiro_experimental afiro --experimental-environment)
if(TEST netlib_afiro_experimental)
    set_tests_properties(
        netlib_afiro_experimental
        PROPERTIES
            ENVIRONMENT
            "LP_SIMPLEX_DISABLE_PRESOLVE=1;LP_SIMPLEX_DISABLE_SCALING=1;LP_SIMPLEX_PROFILE=1"
    )
endif()

# Keep numerical regressions individually named for diagnostics.
foreach(
    model IN
    ITEMS degen2
          25fv47
          agg2
          bnl1
          etamacro
          grow22
          maros
          perold
          pilot.ja
          pilot.we
          pilot87
          qap8
          scfxm1
          share1b
          tuff)
    add_netlib_test(netlib_${model}_dual_revised ${model} --algorithm
                    dual-revised)
endforeach()
