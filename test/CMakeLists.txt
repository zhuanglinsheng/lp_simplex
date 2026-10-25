# Numerical test sources call libm directly; shared-library private linkage
# does not satisfy those references on Linux linkers.
if(UNIX AND NOT APPLE)
    link_libraries(m)
endif()

add_library(test_problem_header OBJECT test_problem_header.c)
target_include_directories(
    test_problem_header
    PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/runtime
            ${PROJECT_SOURCE_DIR}/src/matrix)

add_executable(test_simple_examples test_simple_examples.c)
target_compile_definitions(
    test_simple_examples
    PRIVATE SIMPLE_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/data/simple_examples")
target_link_libraries(test_simple_examples PRIVATE lp_simplex)
add_test(NAME test_simple_examples COMMAND test_simple_examples)

add_executable(test_degeneracy test_degeneracy.c)
target_include_directories(test_degeneracy
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_degeneracy PRIVATE lp_simplex)
add_test(NAME test_degeneracy COMMAND test_degeneracy)

add_executable(test_dual_pricing test_dual_pricing.c)
target_include_directories(test_dual_pricing
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_dual_pricing PRIVATE lp_simplex)
add_test(NAME test_dual_pricing COMMAND test_dual_pricing)

add_executable(test_sparse_vector test_sparse_vector.c)
target_include_directories(test_sparse_vector
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_sparse_vector PRIVATE lp_simplex)
add_test(NAME test_sparse_vector COMMAND test_sparse_vector)

add_executable(test_pan_pdhg_kernels test_pan_pdhg_kernels.c)
target_include_directories(test_pan_pdhg_kernels
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_pdhg_kernels PRIVATE lp_simplex)
add_test(NAME test_pan_pdhg_kernels COMMAND test_pan_pdhg_kernels)

add_executable(test_dual_seed test_dual_seed.c)
target_include_directories(test_dual_seed
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_dual_seed PRIVATE lp_simplex)
add_test(NAME test_dual_seed COMMAND test_dual_seed)
add_test(NAME test_dual_seed_unscaled COMMAND test_dual_seed)
set_tests_properties(test_dual_seed_unscaled
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_DISABLE_SCALING=1")

add_executable(test_cover test_cover.c)
target_include_directories(test_cover
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_cover PRIVATE lp_simplex)
add_test(NAME test_cover COMMAND test_cover)
add_test(NAME test_cover_unscaled COMMAND test_cover)
set_tests_properties(test_cover_unscaled
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_DISABLE_SCALING=1")

add_executable(test_pan test_pan.c)
add_executable(test_pan_bounded test_pan_bounded.c)
target_include_directories(test_pan_bounded
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_bounded PRIVATE lp_simplex)
add_test(NAME test_pan_bounded COMMAND test_pan_bounded)
set_tests_properties(test_pan_bounded PROPERTIES ENVIRONMENT
                                                 "LP_SIMPLEX_PAN_CONE=1")
add_test(NAME test_pan_bounded_implicit COMMAND test_pan_bounded)
set_tests_properties(test_pan_bounded_implicit
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=implicit")
add_test(NAME test_pan_native_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_native_known_optima
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=1")
add_test(NAME test_pan_native_warmstart COMMAND test_pan_warmstart)
set_tests_properties(
    test_pan_native_warmstart
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=1")
add_test(NAME test_pan_native_models COMMAND test_pan)
set_tests_properties(
    test_pan_native_models
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=1")
add_test(NAME test_pan_pdhg_bounded COMMAND test_pan_bounded)
set_tests_properties(
    test_pan_pdhg_bounded
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_NATIVE_PDHG=10000;LP_SIMPLEX_PAN_PDHG_RESTORE=1"
)
add_test(NAME test_pan_pdhg_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_pdhg_known_optima
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_NATIVE_PDHG=10000;LP_SIMPLEX_PAN_PDHG_RESTORE=1"
)
add_test(NAME test_pan_pdhg_models COMMAND test_pan)
set_tests_properties(
    test_pan_pdhg_models
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_NATIVE_PDHG=10000;LP_SIMPLEX_PAN_PDHG_RESTORE=1"
)
add_test(NAME test_pan_native_adaptive COMMAND test_pan_bounded)
set_tests_properties(
    test_pan_native_adaptive
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_NATIVE_ADAPTIVE=1;LP_SIMPLEX_PAN_NATIVE_PREDICT=1"
)
target_link_libraries(test_pan PRIVATE lp_simplex)
add_test(NAME test_pan COMMAND test_pan)

add_executable(test_pan_dynamic_sparse_qr test_pan_dynamic_sparse_qr.c)
target_include_directories(test_pan_dynamic_sparse_qr
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_dynamic_sparse_qr PRIVATE lp_simplex)
add_test(NAME test_pan_dynamic_sparse_qr COMMAND test_pan_dynamic_sparse_qr)

if(LP_SIMPLEX_USE_SPQR
   AND SPQR_INCLUDE_DIR
   AND SPQR_LIBRARY)
    add_executable(test_pan_multifrontal_qr test_pan_multifrontal_qr.c)
    target_include_directories(test_pan_multifrontal_qr
                               PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
    target_include_directories(test_pan_multifrontal_qr SYSTEM
                               PRIVATE ${SPQR_INCLUDE_DIR})
    target_link_libraries(
        test_pan_multifrontal_qr PRIVATE lp_simplex ${CHOLMOD_LIBRARY}
                                         ${SUITESPARSE_CONFIG_LIBRARY})
    add_test(NAME test_pan_multifrontal_qr COMMAND test_pan_multifrontal_qr)
endif()

add_executable(test_presolve test_presolve.c)
target_include_directories(test_presolve
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_presolve PRIVATE lp_simplex)
add_test(NAME test_presolve COMMAND test_presolve)

add_executable(test_presolve_dominance test_presolve_dominance.c)
target_include_directories(test_presolve_dominance
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_presolve_dominance PRIVATE lp_simplex)
if(NOT MSVC)
    target_link_libraries(test_presolve_dominance PRIVATE m)
endif()
add_test(NAME test_presolve_dominance COMMAND test_presolve_dominance)

add_executable(test_scaling test_scaling.c)
target_include_directories(test_scaling
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_scaling PRIVATE lp_simplex)
add_test(NAME test_scaling COMMAND test_scaling)

add_executable(test_basis_update test_basis_update.c)
target_include_directories(test_basis_update
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_basis_update PRIVATE lp_simplex)
add_test(NAME test_basis_update COMMAND test_basis_update)

add_executable(test_regressions test_regressions.c)
target_compile_definitions(
    test_regressions
    PRIVATE SIMPLE_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/data/simple_examples")
target_link_libraries(test_regressions PRIVATE lp_simplex)
add_test(NAME test_regressions COMMAND test_regressions)

add_executable(test_primal_pricing test_primal_pricing.c)
target_include_directories(test_primal_pricing
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_primal_pricing PRIVATE lp_simplex)
add_test(NAME test_primal_pricing COMMAND test_primal_pricing)
add_test(NAME test_regressions_devex COMMAND test_regressions)
set_tests_properties(test_regressions_devex
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_CLEANUP_DEVEX=1")

add_executable(test_pan_phase1 test_pan_phase1.c)
target_include_directories(test_pan_phase1
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_phase1 PRIVATE lp_simplex)
add_test(NAME test_pan_phase1 COMMAND test_pan_phase1)
set_tests_properties(test_pan_phase1
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_SINGLETON_CRASH=1")

add_test(NAME test_regressions_initial_cost_perturb COMMAND test_regressions)
set_tests_properties(test_regressions_initial_cost_perturb
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_INITIAL_COST_PERTURB=1")
add_test(NAME test_regressions_stable_klu COMMAND test_regressions)
set_tests_properties(test_regressions_stable_klu
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_STABLE_KLU=1")

add_executable(test_recovery test_recovery.c)
target_include_directories(test_recovery
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_recovery PRIVATE lp_simplex)
add_test(NAME test_recovery COMMAND test_recovery)

add_executable(test_dual_recovery_restart test_dual_recovery_restart.c)
target_include_directories(test_dual_recovery_restart
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_dual_recovery_restart PRIVATE lp_simplex)
add_test(NAME test_dual_recovery_restart COMMAND test_dual_recovery_restart)
add_test(NAME test_regressions_bound_recovery COMMAND test_regressions)
set_tests_properties(
    test_regressions_bound_recovery
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_FORCE_BOUND_PERTURB=1;LP_SIMPLEX_INITIAL_COST_PERTURB=1")
add_test(NAME test_regressions_auto_bound_recovery COMMAND test_regressions)
set_tests_properties(test_regressions_auto_bound_recovery
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PRIMAL_BOUND_PERTURB=1")
add_test(NAME test_regressions_no_bound_perturb COMMAND test_regressions)
set_tests_properties(
    test_regressions_no_bound_perturb
    PROPERTIES ENVIRONMENT "LP_SIMPLEX_DISABLE_BOUND_PERTURB=1")

add_executable(test_pan_warmstart test_pan_warmstart.c)
target_include_directories(test_pan_warmstart
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_warmstart PRIVATE lp_simplex)
add_test(NAME test_pan_warmstart COMMAND test_pan_warmstart)

if(LP_SIMPLEX_USE_SPQR
   AND SPQR_INCLUDE_DIR
   AND SPQR_LIBRARY)
    add_executable(test_pan_qr_stress test_pan_qr_stress.c)
    target_include_directories(test_pan_qr_stress
                               PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
    target_include_directories(test_pan_qr_stress SYSTEM
                               PRIVATE ${SPQR_INCLUDE_DIR})
    target_link_libraries(
        test_pan_qr_stress PRIVATE lp_simplex ${CHOLMOD_LIBRARY}
                                   ${SUITESPARSE_CONFIG_LIBRARY})
    add_test(NAME test_pan_qr_stress COMMAND test_pan_qr_stress)
endif()

if(LP_SIMPLEX_USE_SPQR
   AND SPQR_INCLUDE_DIR
   AND SPQR_LIBRARY)
    add_executable(test_pan_rank_repair test_pan_rank_repair.c)
    target_include_directories(test_pan_rank_repair
                               PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
    target_link_libraries(test_pan_rank_repair PRIVATE lp_simplex)
    add_test(NAME test_pan_rank_repair COMMAND test_pan_rank_repair)
endif()

add_executable(test_pan_basis_stress test_pan_basis_stress.c)
target_include_directories(test_pan_basis_stress
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_basis_stress PRIVATE lp_simplex)
add_test(NAME test_pan_basis_stress COMMAND test_pan_basis_stress)
if(LP_SIMPLEX_USE_SPQR
   AND SPQR_INCLUDE_DIR
   AND SPQR_LIBRARY)
    add_test(NAME test_pan_large_sparse COMMAND test_pan_basis_stress --large)
    set_tests_properties(test_pan_large_sparse PROPERTIES TIMEOUT 90 LABELS
                                                          "large-sparse")
endif()

add_executable(test_pan_phase1_accuracy test_pan_phase1_accuracy.c)
target_include_directories(test_pan_phase1_accuracy
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_phase1_accuracy PRIVATE lp_simplex)
add_test(NAME test_pan_phase1_accuracy COMMAND test_pan_phase1_accuracy)

add_executable(test_pan_known_optima test_pan_known_optima.c)
target_link_libraries(test_pan_known_optima PRIVATE lp_simplex)
add_test(NAME test_pan_known_optima COMMAND test_pan_known_optima)

add_executable(test_pan_mixed_nnls test_pan_mixed_nnls.c)
target_include_directories(test_pan_mixed_nnls
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_mixed_nnls PRIVATE lp_simplex)
add_test(NAME test_pan_mixed_nnls COMMAND test_pan_mixed_nnls)
add_executable(test_pan_cone_oracle test_pan_cone_oracle.c)
add_executable(test_pan_trial test_pan_trial.c)
target_include_directories(test_pan_trial
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_trial PRIVATE lp_simplex)
add_test(NAME test_pan_trial COMMAND test_pan_trial)
add_executable(test_pan_recovery test_pan_recovery.c)
target_include_directories(test_pan_recovery
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_recovery PRIVATE lp_simplex)
add_test(NAME test_pan_recovery COMMAND test_pan_recovery)
add_test(NAME test_pan_recovery_implicit COMMAND test_pan_recovery)
set_tests_properties(test_pan_recovery_implicit
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=implicit")
add_test(NAME test_pan_recovery_centered COMMAND test_pan_recovery)
set_tests_properties(
    test_pan_recovery_centered
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_DUAL_CENTER=1;LP_SIMPLEX_PAN_DUAL_POLISH=128;LP_SIMPLEX_PAN_CENTER_EPOCH=2"
)
add_test(NAME test_pan_recovery_predictor COMMAND test_pan_recovery)
set_tests_properties(
    test_pan_recovery_predictor
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_FIRST_ORDER=1;LP_SIMPLEX_PAN_DUAL_POLISH=1;LP_SIMPLEX_PAN_FACE_GUARD=1"
)
target_include_directories(test_pan_cone_oracle
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_cone_oracle PRIVATE lp_simplex)
add_test(NAME test_pan_cone_oracle COMMAND test_pan_cone_oracle)
add_test(NAME test_pan_cone_known_optima COMMAND test_pan_known_optima)
set_tests_properties(test_pan_cone_known_optima
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=1")
add_test(NAME test_pan_cone_warmstart COMMAND test_pan_warmstart)
set_tests_properties(test_pan_cone_warmstart PROPERTIES ENVIRONMENT
                                                        "LP_SIMPLEX_PAN_CONE=1")
add_test(NAME test_pan_cone_models COMMAND test_pan)
set_tests_properties(test_pan_cone_models PROPERTIES ENVIRONMENT
                                                     "LP_SIMPLEX_PAN_CONE=1")
add_test(NAME test_pan_implicit_known_optima COMMAND test_pan_known_optima)
set_tests_properties(test_pan_implicit_known_optima
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=implicit")
add_test(NAME test_pan_implicit_warmstart COMMAND test_pan_warmstart)
set_tests_properties(test_pan_implicit_warmstart
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=implicit")
add_test(NAME test_pan_implicit_models COMMAND test_pan)
set_tests_properties(test_pan_implicit_models
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=implicit")
add_test(NAME test_pan_centered_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_centered_known_optima
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_DUAL_CENTER=1")
add_test(NAME test_pan_centered_warmstart COMMAND test_pan_warmstart)
set_tests_properties(
    test_pan_centered_warmstart
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_DUAL_CENTER=1")
add_test(NAME test_pan_centered_models COMMAND test_pan)
set_tests_properties(
    test_pan_centered_models
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_DUAL_CENTER=1")
add_test(NAME test_pan_batch_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_batch_known_optima
    PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_CONE_BATCH=1")

add_test(NAME test_pan_block_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_block_known_optima
    PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_NNLS_BLOCK=16")
add_test(NAME test_pan_face_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_face_known_optima
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_NNLS_BLOCK=16;LP_SIMPLEX_PAN_FACE_GUARD=1"
)
add_test(NAME test_pan_face_implicit_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_face_implicit_known_optima
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_FACE_GUARD=1")
add_test(NAME test_pan_first_order_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_first_order_known_optima
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_NNLS_BLOCK=16;LP_SIMPLEX_PAN_FACE_GUARD=1;LP_SIMPLEX_PAN_FIRST_ORDER=2048"
)
add_test(NAME test_pan_first_order_warmstart COMMAND test_pan_warmstart)
set_tests_properties(
    test_pan_first_order_warmstart
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_NNLS_BLOCK=16;LP_SIMPLEX_PAN_FACE_GUARD=1;LP_SIMPLEX_PAN_FIRST_ORDER=2048"
)
add_test(NAME test_pan_dual_polish_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_dual_polish_known_optima
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_DUAL_POLISH=2048")
add_test(NAME test_pan_epoch_known_optima COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_epoch_known_optima
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_DUAL_CENTER=1;LP_SIMPLEX_PAN_DUAL_POLISH=2048;LP_SIMPLEX_PAN_CENTER_EPOCH=2;LP_SIMPLEX_PAN_FACE_GUARD=1;LP_SIMPLEX_PAN_NNLS_BLOCK=16;LP_SIMPLEX_PAN_FIRST_ORDER=2048"
)
add_test(NAME test_pan_epoch_warmstart COMMAND test_pan_warmstart)
set_tests_properties(
    test_pan_epoch_warmstart
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_DUAL_CENTER=1;LP_SIMPLEX_PAN_DUAL_POLISH=2048;LP_SIMPLEX_PAN_CENTER_EPOCH=2;LP_SIMPLEX_PAN_FACE_GUARD=1;LP_SIMPLEX_PAN_NNLS_BLOCK=16;LP_SIMPLEX_PAN_FIRST_ORDER=2048"
)
add_test(NAME test_pan_implicit_updates_known_optima
         COMMAND test_pan_known_optima)
set_tests_properties(
    test_pan_implicit_updates_known_optima
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_IMPLICIT_UPDATES=1")
add_test(NAME test_pan_implicit_updates_warmstart COMMAND test_pan_warmstart)
set_tests_properties(
    test_pan_implicit_updates_warmstart
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_IMPLICIT_UPDATES=1")
add_test(NAME test_pan_implicit_updates_models COMMAND test_pan)
set_tests_properties(
    test_pan_implicit_updates_models
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_IMPLICIT_UPDATES=1")
if(LP_SIMPLEX_USE_SPQR
   AND SPQR_INCLUDE_DIR
   AND SPQR_LIBRARY)
    add_executable(test_pan_face test_pan_face.c)
    target_include_directories(test_pan_face
                               PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
    target_link_libraries(test_pan_face PRIVATE lp_simplex)
    add_test(NAME test_pan_face COMMAND test_pan_face)
    set_tests_properties(
        test_pan_face
        PROPERTIES ENVIRONMENT
                   "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_FACE_GUARD=1")
    add_test(NAME test_pan_face_implicit COMMAND test_pan_face)
    set_tests_properties(
        test_pan_face_implicit
        PROPERTIES ENVIRONMENT
                   "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_FACE_GUARD=1")
    add_test(NAME test_pan_face_predictor COMMAND test_pan_face)
    set_tests_properties(
        test_pan_face_predictor
        PROPERTIES
            ENVIRONMENT
            "LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_FACE_GUARD=1;LP_SIMPLEX_PAN_FIRST_ORDER=1;LP_SIMPLEX_PAN_DUAL_POLISH=1"
    )
endif()

add_executable(test_auto test_auto.c)
target_include_directories(test_auto PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto PRIVATE lp_simplex)
add_test(NAME test_auto COMMAND test_auto)
add_executable(test_auto_known_optima test_pan_known_optima.c)
target_compile_definitions(test_auto_known_optima PRIVATE TEST_AUTO=1)
target_link_libraries(test_auto_known_optima PRIVATE lp_simplex)
add_test(NAME test_auto_known_optima COMMAND test_auto_known_optima)
add_executable(test_auto_flow test_auto_flow.c)
target_include_directories(test_auto_flow
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_flow PRIVATE lp_simplex)
add_test(NAME test_auto_flow COMMAND test_auto_flow)
add_test(NAME test_auto_settings COMMAND test_auto)
set_tests_properties(
    test_auto_settings
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_CONE=1;LP_SIMPLEX_PAN_NATIVE_PDHG=0;LP_SIMPLEX_PAN_DUAL_POLISH=0;LP_SIMPLEX_PAN_NATIVE_PREDICT=1;LP_SIMPLEX_PAN_IMPLICIT_BULK=1;LP_SIMPLEX_PAN_FULL_REFINEMENT=1;LP_SIMPLEX_PAN_DISABLE_SINGLETON_CRASH=1;LP_SIMPLEX_PAN_BLOCK_CRASH=1;LP_SIMPLEX_PAN_NNLS_BLOCK=32;LP_SIMPLEX_PAN_DUAL_CENTER=1;LP_SIMPLEX_PAN_NEWTON_QR=1"
)
# These variants deliberately share lp_simplex_regression.mps in the build
# directory. Parallel writers can change fixture dimensions during a read.
set_tests_properties(
    test_regressions
    test_regressions_devex
    test_regressions_initial_cost_perturb
    test_regressions_stable_klu
    test_regressions_bound_recovery
    test_regressions_auto_bound_recovery
    test_regressions_no_bound_perturb
    PROPERTIES RESOURCE_LOCK regression_mps)
add_executable(test_result_contract test_result_contract.c)
target_include_directories(test_result_contract
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_result_contract PRIVATE lp_simplex)
add_test(NAME test_result_contract COMMAND test_result_contract)
add_test(NAME test_result_contract_native COMMAND test_result_contract)
set_tests_properties(
    test_result_contract_native
    PROPERTIES ENVIRONMENT
               "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=implicit")
add_test(NAME test_result_contract_cone COMMAND test_result_contract)
set_tests_properties(test_result_contract_cone
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=implicit")
add_test(NAME test_result_contract_unscaled COMMAND test_result_contract)
set_tests_properties(test_result_contract_unscaled
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_DISABLE_SCALING=1")

add_executable(test_auto_policy test_auto_policy.c)
target_include_directories(test_auto_policy
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_policy PRIVATE lp_simplex)
add_test(NAME test_auto_policy COMMAND test_auto_policy)
add_test(NAME test_auto_policy_settings COMMAND test_auto_policy)
set_tests_properties(
    test_auto_policy_settings
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_DUAL_START=1;LP_SIMPLEX_PAN_COMPLETE_START=1;LP_SIMPLEX_PAN_COMPATIBLE_PRICING=1;LP_SIMPLEX_PAN_STABILIZE=1;LP_SIMPLEX_PAN_FIRST_ORDER=1000"
)

# Optional manual original-model certificate export for numerical audits.
add_executable(pan_native_certificate EXCLUDE_FROM_ALL pan_native_certificate.c)
target_include_directories(pan_native_certificate
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(pan_native_certificate PRIVATE lp_simplex)

add_executable(test_progress test_progress.c)
target_include_directories(test_progress
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_progress PRIVATE lp_simplex)
add_test(NAME test_progress COMMAND test_progress)
add_test(NAME test_progress_native COMMAND test_progress)
set_tests_properties(
    test_progress_native
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_NATIVE_PDHG=10000;LP_SIMPLEX_PAN_CONE=implicit"
)

add_executable(test_infeasibility_certificate test_infeasibility_certificate.c)
target_include_directories(test_infeasibility_certificate
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_infeasibility_certificate PRIVATE lp_simplex)
if(LP_SIMPLEX_USE_GMP
   AND GMP_INCLUDE_DIR
   AND GMP_LIBRARY)
    target_compile_definitions(test_infeasibility_certificate
                               PRIVATE TEST_HAVE_GMP)
endif()
add_test(NAME test_infeasibility_certificate
         COMMAND test_infeasibility_certificate)

add_executable(test_numerical_policy test_numerical_policy.c)
target_include_directories(test_numerical_policy
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_numerical_policy PRIVATE lp_simplex)
add_test(NAME test_numerical_policy COMMAND test_numerical_policy)
add_test(NAME test_numerical_policy_pdhg COMMAND test_numerical_policy)
set_tests_properties(
    test_numerical_policy_pdhg
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_NATIVE_PDHG=10000;LP_SIMPLEX_PAN_PDHG_RESTORE=1"
)

add_executable(test_stage_report test_stage_report.c)
target_include_directories(test_stage_report
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_stage_report PRIVATE lp_simplex)
add_test(NAME test_stage_report COMMAND test_stage_report)

add_executable(test_execution test_execution.c)
target_include_directories(test_execution
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_execution PRIVATE lp_simplex)
add_test(NAME test_execution COMMAND test_execution)

add_executable(test_checkpoint test_checkpoint.c)
target_include_directories(test_checkpoint
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_checkpoint PRIVATE lp_simplex)
add_test(NAME test_checkpoint COMMAND test_checkpoint)

add_executable(test_retry test_retry.c)
target_include_directories(test_retry
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_retry PRIVATE lp_simplex)
add_test(NAME test_retry COMMAND test_retry)
add_test(NAME test_retry_unscaled COMMAND test_retry)
set_tests_properties(test_retry_unscaled
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_DISABLE_SCALING=1")

add_test(NAME test_retry_overrides COMMAND test_retry)
set_tests_properties(test_retry_overrides
                     PROPERTIES ENVIRONMENT "LP_SIMPLEX_PAN_CONE=implicit")

add_executable(test_auto_state test_auto_state.c)
target_include_directories(test_auto_state
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_state PRIVATE lp_simplex)
add_test(NAME test_auto_state COMMAND test_auto_state)

add_executable(test_bda_config test_bda_config.c)
target_include_directories(test_bda_config
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_bda_config PRIVATE lp_simplex)
add_test(NAME test_bda_config COMMAND test_bda_config)
add_test(NAME test_bda_config_settings COMMAND test_bda_config)
set_tests_properties(
    test_bda_config_settings
    PROPERTIES
        ENVIRONMENT
        "LP_SIMPLEX_PAN_NATIVE_BOUNDS=1;LP_SIMPLEX_PAN_CONE=implicit;LP_SIMPLEX_PAN_NATIVE_PDHG=0;LP_SIMPLEX_PAN_FIRST_ORDER=999;LP_SIMPLEX_PAN_DUAL_POLISH=0;LP_SIMPLEX_PAN_DUAL_CENTER=1;LP_SIMPLEX_PAN_CENTER_EPOCH=3;LP_SIMPLEX_PAN_CONE_BATCH=1;LP_SIMPLEX_PAN_FACE_GUARD=1;LP_SIMPLEX_PAN_NATIVE_PREDICT=1;LP_SIMPLEX_PAN_NNLS_BLOCK=32;LP_SIMPLEX_PAN_IMPLICIT_BULK=1;LP_SIMPLEX_PAN_FULL_REFINEMENT=1;LP_SIMPLEX_PAN_DISABLE_SINGLETON_CRASH=1"
)

add_executable(test_pan_minimum_norm test_pan_minimum_norm.c)
target_include_directories(test_pan_minimum_norm
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_minimum_norm PRIVATE lp_simplex)
if(NOT
   (LP_SIMPLEX_USE_SPQR
    AND SPQR_INCLUDE_DIR
    AND SPQR_LIBRARY
    AND CHOLMOD_LIBRARY
    AND SUITESPARSE_CONFIG_LIBRARY))
    target_compile_definitions(test_pan_minimum_norm
                               PRIVATE TEST_DENSE_BACKEND=1)
endif()
add_test(NAME test_pan_minimum_norm COMMAND test_pan_minimum_norm)

add_executable(test_pan_services test_pan_services.c)
target_include_directories(test_pan_services
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_services PRIVATE lp_simplex)
add_test(NAME test_pan_services COMMAND test_pan_services)

add_executable(test_auto_contract test_auto_contract.c)
target_include_directories(test_auto_contract
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_contract PRIVATE lp_simplex)
add_test(NAME test_auto_contract COMMAND test_auto_contract)

add_executable(test_column_handoff test_column_handoff.c)
target_include_directories(test_column_handoff
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_column_handoff PRIVATE lp_simplex)
add_test(NAME test_column_handoff COMMAND test_column_handoff)

add_executable(test_environment test_environment.c)
target_include_directories(test_environment
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_environment PRIVATE lp_simplex)
add_test(NAME test_environment COMMAND test_environment)

add_executable(test_basis_retry test_basis_retry.c)
target_include_directories(test_basis_retry
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_basis_retry PRIVATE lp_simplex)
add_test(NAME test_basis_retry COMMAND test_basis_retry)

add_executable(test_auto_output test_auto_output.c)
target_include_directories(test_auto_output
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_output PRIVATE lp_simplex)
add_test(NAME test_auto_output COMMAND test_auto_output)

add_executable(test_progress_origin test_progress_origin.c)
target_include_directories(test_progress_origin
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_progress_origin PRIVATE lp_simplex)
add_test(NAME test_progress_origin COMMAND test_progress_origin)

add_executable(test_postsolve test_postsolve.c)
target_include_directories(test_postsolve
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_postsolve PRIVATE lp_simplex)
add_test(NAME test_postsolve COMMAND test_postsolve)

add_executable(test_postsolve_classification test_postsolve_classification.c)
target_include_directories(test_postsolve_classification
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_postsolve_classification PRIVATE lp_simplex)
add_test(NAME test_postsolve_classification
         COMMAND test_postsolve_classification)

add_executable(test_auto_initializer test_auto_initializer.c)
target_include_directories(test_auto_initializer
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_initializer PRIVATE lp_simplex)
add_test(NAME test_auto_initializer COMMAND test_auto_initializer)

add_executable(test_basis_origin test_basis_origin.c)
target_include_directories(test_basis_origin
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_basis_origin PRIVATE lp_simplex)
add_test(NAME test_basis_origin COMMAND test_basis_origin)

add_executable(test_checkpoint_fault test_checkpoint_fault.c)
target_include_directories(test_checkpoint_fault
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_compile_definitions(
    test_checkpoint_fault
    PRIVATE $<TARGET_PROPERTY:lp_simplex,COMPILE_DEFINITIONS>)
target_link_libraries(test_checkpoint_fault
                      PRIVATE lp_simplex ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
add_test(NAME test_checkpoint_fault COMMAND test_checkpoint_fault)

add_executable(test_retry_origin test_retry_origin.c)
target_include_directories(test_retry_origin
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_retry_origin PRIVATE lp_simplex)
add_test(NAME test_retry_origin COMMAND test_retry_origin)

add_executable(test_auto_capabilities test_auto_capabilities.c)
target_include_directories(test_auto_capabilities
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_capabilities PRIVATE lp_simplex)
add_test(NAME test_auto_capabilities COMMAND test_auto_capabilities)

# Installed-package consumers build in an independent project after relocation.
# Python tooling is optional for C-only users; CI installs its test dependency.
include(CheckLanguage)
check_language(CXX)
if(CMAKE_CXX_COMPILER)
    set(LP_SIMPLEX_TEST_CXX ON)
else()
    set(LP_SIMPLEX_TEST_CXX OFF)
    message(
        STATUS "C++ compiler unavailable; installed C consumer remains enabled")
endif()
if(CMAKE_C_FLAGS_DEBUG MATCHES "fsanitize=address")
    set(LP_SIMPLEX_DEBUG_ASAN ON)
else()
    set(LP_SIMPLEX_DEBUG_ASAN OFF)
endif()
if(NOT CMAKE_CROSSCOMPILING)
    add_test(
        NAME test_install_package
        COMMAND
            ${CMAKE_COMMAND} "-DBUILD_DIR=${PROJECT_BINARY_DIR}"
            "-DWORK_DIR=${PROJECT_BINARY_DIR}/test/installed-package"
            "-DPROJECT_SOURCE=${PROJECT_SOURCE_DIR}"
            "-DCONSUMER_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}/package_consumer"
            "-DCTEST_COMMAND=${CMAKE_CTEST_COMMAND}" "-DCONFIG=$<CONFIG>"
            "-DGENERATOR=${CMAKE_GENERATOR}"
            "-DGENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}"
            "-DGENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}"
            "-DC_COMPILER=${CMAKE_C_COMPILER}"
            "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DTEST_CXX=${LP_SIMPLEX_TEST_CXX}"
            "-DDEBUG_ASAN=${LP_SIMPLEX_DEBUG_ASAN}"
            "-DDATAROOT=${CMAKE_INSTALL_DATAROOTDIR}"
            "-DBINDIR=${CMAKE_INSTALL_BINDIR}"
            "-DLIBDIR=${CMAKE_INSTALL_LIBDIR}" -P
            ${CMAKE_CURRENT_SOURCE_DIR}/test_install_package.cmake)
    set_tests_properties(test_install_package PROPERTIES LABELS "package")
endif()
find_package(
    Python3
    COMPONENTS Interpreter
    QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME test_architecture
             COMMAND ${Python3_EXECUTABLE}
                     ${PROJECT_SOURCE_DIR}/tools/check_architecture.py)
    add_test(NAME test_prepare_data
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/test_prepare_data.py)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy"
        RESULT_VARIABLE LP_SIMPLEX_TEST_PYTHON_NUMERICS
        OUTPUT_QUIET ERROR_QUIET)
    if(LP_SIMPLEX_TEST_PYTHON_NUMERICS EQUAL 0)
        add_test(NAME test_benchmark_tools
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/test_benchmark_solvers.py)
    else()
        message(
            STATUS
                "NumPy/SciPy unavailable; benchmark-tool Python regression not registered"
        )
    endif()
endif()

add_executable(test_tableau_budget test_tableau_budget.c)
target_include_directories(test_tableau_budget
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_tableau_budget PRIVATE lp_simplex)
add_test(NAME test_tableau_budget COMMAND test_tableau_budget)

add_executable(test_tableau_report test_tableau_report.c)
target_include_directories(test_tableau_report
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_tableau_report PRIVATE lp_simplex)
add_test(NAME test_tableau_report COMMAND test_tableau_report)

add_executable(test_tableau_model test_tableau_model.c)
target_include_directories(test_tableau_model
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_tableau_model PRIVATE lp_simplex)
add_test(NAME test_tableau_model COMMAND test_tableau_model)

add_executable(test_tableau_stability test_tableau_stability.c)
target_include_directories(test_tableau_stability
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_tableau_stability PRIVATE lp_simplex)
add_test(NAME test_tableau_stability COMMAND test_tableau_stability)

if(Python3_Interpreter_FOUND AND NOT CMAKE_CROSSCOMPILING)
    # Keep this numerical regression available when benchmark targets are off.
    add_executable(test_tableau_netlib_driver
                   ${PROJECT_SOURCE_DIR}/data/lp_solve_mps.c)
    target_link_libraries(test_tableau_netlib_driver PRIVATE lp_simplex)
    add_test(
        NAME test_tableau_netlib
        COMMAND
            ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_tableau_netlib.py
            $<TARGET_FILE:test_tableau_netlib_driver>)
    set_tests_properties(test_tableau_netlib PROPERTIES TIMEOUT 90)
endif()

add_executable(test_tableau_storage test_tableau_storage.c)
target_include_directories(test_tableau_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_tableau_storage PRIVATE lp_simplex)
add_test(NAME test_tableau_storage COMMAND test_tableau_storage)

add_executable(test_singleton_classification test_singleton_classification.c)
target_include_directories(test_singleton_classification
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_singleton_classification PRIVATE lp_simplex)
add_test(NAME test_singleton_classification
         COMMAND test_singleton_classification)

add_executable(test_feasibility_retry test_feasibility_retry.c)
target_include_directories(test_feasibility_retry
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_feasibility_retry PRIVATE lp_simplex)
add_test(NAME test_feasibility_retry COMMAND test_feasibility_retry)

add_executable(test_model_storage test_model_storage.c)
target_include_directories(test_model_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_model_storage PRIVATE lp_simplex)
add_test(NAME test_model_storage COMMAND test_model_storage)

add_executable(test_pan_storage test_pan_storage.c)
target_include_directories(test_pan_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_storage PRIVATE lp_simplex)
add_test(NAME test_pan_storage COMMAND test_pan_storage)

add_executable(test_pan_qr_storage test_pan_qr_storage.c)
target_include_directories(test_pan_qr_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_qr_storage PRIVATE lp_simplex)
add_test(NAME test_pan_qr_storage COMMAND test_pan_qr_storage)

add_executable(test_dual_storage test_dual_storage.c)
target_include_directories(test_dual_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_dual_storage PRIVATE lp_simplex)
add_test(NAME test_dual_storage COMMAND test_dual_storage)
if(LP_SIMPLEX_USE_KLU
   AND KLU_INCLUDE_DIR
   AND KLU_LIBRARY)
    add_executable(test_dual_storage_klu test_dual_storage.c)
    target_compile_definitions(test_dual_storage_klu
                               PRIVATE LP_SIMPLEX_HAVE_KLU)
    target_include_directories(test_dual_storage_klu
                               PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
    target_include_directories(test_dual_storage_klu SYSTEM
                               PRIVATE ${KLU_INCLUDE_DIR})
    target_link_libraries(test_dual_storage_klu PRIVATE lp_simplex
                                                        ${KLU_LIBRARY})
    add_test(NAME test_dual_storage_klu COMMAND test_dual_storage_klu)
endif()

add_executable(test_sparse_lu_storage test_sparse_lu_storage.c)
target_include_directories(test_sparse_lu_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_sparse_lu_storage PRIVATE lp_simplex)
add_test(NAME test_sparse_lu_storage COMMAND test_sparse_lu_storage)

add_executable(test_problem_storage test_problem_storage.c)
target_include_directories(test_problem_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_problem_storage PRIVATE lp_simplex)
add_test(NAME test_problem_storage COMMAND test_problem_storage)

add_executable(test_presolve_storage test_presolve_storage.c)
target_include_directories(test_presolve_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_presolve_storage PRIVATE lp_simplex)
add_test(NAME test_presolve_storage COMMAND test_presolve_storage)

add_executable(test_presolve_record test_presolve_record.c)
target_include_directories(test_presolve_record
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_presolve_record PRIVATE lp_simplex)
add_test(NAME test_presolve_record COMMAND test_presolve_record)

add_executable(test_dominance_storage test_dominance_storage.c)
target_include_directories(test_dominance_storage
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_dominance_storage PRIVATE lp_simplex)
add_test(NAME test_dominance_storage COMMAND test_dominance_storage)

add_executable(test_dominance_storage_small test_dominance_storage.c)
target_include_directories(test_dominance_storage_small
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_compile_definitions(test_dominance_storage_small
                           PRIVATE TEST_DOMINANCE_BYTES=256)
target_link_libraries(test_dominance_storage_small PRIVATE lp_simplex)
add_test(NAME test_dominance_storage_small COMMAND test_dominance_storage_small)

add_executable(test_pdhg_restore_trigger test_pdhg_restore_trigger.c)
target_include_directories(test_pdhg_restore_trigger
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pdhg_restore_trigger PRIVATE lp_simplex)
add_test(NAME test_pdhg_restore_trigger COMMAND test_pdhg_restore_trigger)

add_executable(test_pan_face_recovery test_pan_face_recovery.c)
target_include_directories(test_pan_face_recovery
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_face_recovery PRIVATE lp_simplex)
add_test(NAME test_pan_face_recovery COMMAND test_pan_face_recovery)

add_executable(test_pan_face_service test_pan_face_service.c)
target_include_directories(test_pan_face_service
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_pan_face_service PRIVATE lp_simplex)
add_test(NAME test_pan_face_service COMMAND test_pan_face_service)

add_executable(test_auto_resource test_auto_resource.c)
target_include_directories(test_auto_resource
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})
target_link_libraries(test_auto_resource PRIVATE lp_simplex)
add_test(NAME test_auto_resource COMMAND test_auto_resource)
