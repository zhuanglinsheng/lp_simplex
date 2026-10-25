###############################################################################
# Setting
###############################################################################

cmake_minimum_required(VERSION 3.12)
project(lp_simplex LANGUAGES C)
include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
option(LP_SIMPLEX_BUILD_BENCHMARKS "Build the MPS benchmark driver" ON)
option(LP_SIMPLEX_BUILD_EXAMPLES "Build standalone C API usage examples" ON)
option(LP_SIMPLEX_USE_KLU "Use SuiteSparse KLU for sparse basis factorizations"
       ON)
option(LP_SIMPLEX_USE_GMP
       "Use GMP for exact terminal infeasibility verification" ON)
option(LP_SIMPLEX_USE_SPQR "Use CHOLMOD/SPQR sparse algebra for Pan/BDA" ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_EXTENSIONS OFF)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # Keep C90 syntax while requesting declarations for the C99 math runtime
    # (fma, hypot, isfinite, ...), hidden by glibc in strict C90 mode.
    add_compile_definitions(_ISOC99_SOURCE=1)
endif()
if(MSVC)
    # Force to always compile with W4
    if(CMAKE_C_FLAGS_DEBUG MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS_DEBUG
                             "${CMAKE_C_FLAGS_DEBUG}")
    else()
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /W4")
    endif()
else()
    # Update if necessary
    set(CMAKE_C_FLAGS_DEBUG
        "${CMAKE_C_FLAGS_DEBUG} -Wall -Wno-long-long -pedantic -fsanitize=address"
    )
endif()
##
## Installation setting
##
# The installed library must not retain the original installation prefix.
# Consumer targets supply their own runtime search path for the imported file.

###############################################################################
# BLAS Required
###############################################################################

find_package(BLAS REQUIRED)

if(BLAS_FOUND)
    message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
    set(USE_BLAS TRUE)
else()
    message(WARNING "BLAS not found")
    set(USE_BLAS FALSE)
endif()

###############################################################################
# LAPACK Required
###############################################################################

find_package(LAPACK REQUIRED)

if(LAPACK_FOUND)
    message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
    set(USE_LAPACK TRUE)
else()
    message(WARNING "LAPACK not found")
    set(USE_LAPACK FALSE)
endif()

###############################################################################
# Target
###############################################################################

set(LP_SIMPLEX_PRIVATE_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/src/core
    ${PROJECT_SOURCE_DIR}/src/runtime
    ${PROJECT_SOURCE_DIR}/src/presolve
    ${PROJECT_SOURCE_DIR}/src/dual
    ${PROJECT_SOURCE_DIR}/src/bda
    ${PROJECT_SOURCE_DIR}/src/bda/linalg
    ${PROJECT_SOURCE_DIR}/src/basis
    ${PROJECT_SOURCE_DIR}/src/matrix
    ${PROJECT_SOURCE_DIR}/src/tableau
    ${PROJECT_SOURCE_DIR}/src/common)

set(SOURCES
    src/core/model.c
    src/core/mps.c
    src/core/simplex.c
    src/core/simplex_model_audit.c
    src/core/simplex_original_recovery.c
    src/runtime/simplex_stage.c
    src/runtime/simplex_execution.c
    src/runtime/simplex_solve_context.c
    src/runtime/simplex_run_config.c
    src/runtime/simplex_retry.c
    src/core/simplex_auto.c
    src/core/simplex_auto_policy.c
    src/core/simplex_auto_routes.c
    src/runtime/simplex_problem.c
    src/runtime/simplex_problem_identity.c
    src/runtime/simplex_certificate.c
    src/runtime/simplex_basis_certificate.c
    src/runtime/simplex_progress.c
    src/runtime/simplex_scaling.c
    src/presolve/simplex_presolve.c
    src/presolve/simplex_presolve_dominance.c
    src/presolve/simplex_presolve_activity.c
    src/presolve/simplex_presolve_postsolve.c
    src/presolve/simplex_presolve_profile.c
    src/presolve/simplex_presolve_queue.c
    src/presolve/simplex_presolve_record.c
    src/presolve/simplex_presolve_rules.c
    src/presolve/simplex_presolve_substitution.c
    src/dual/simplex_dual.c
    src/dual/simplex_dual_backend.c
    src/dual/simplex_dual_kernel.c
    src/dual/simplex_dual_checkpoint.c
    src/dual/simplex_column_handoff.c
    src/dual/simplex_dual_retry.c
    src/dual/simplex_cover.c
    src/dual/simplex_dual_bounds.c
    src/dual/simplex_dual_state.c
    src/dual/simplex_dual_feasibility.c
    src/dual/simplex_dual_pricing.c
    src/dual/simplex_dual_primal_pricing.c
    src/dual/simplex_singleton_dual.c
    src/dual/simplex_degeneracy.c
    src/dual/simplex_recovery.c
    src/bda/simplex_pan_service.c
    src/bda/simplex_pan_config.c
    src/bda/simplex_pan.c
    src/bda/simplex_pan_bounded.c
    src/bda/simplex_pan_face_recovery.c
    src/bda/simplex_pan_pdhg.c
    src/bda/simplex_pan_cone.c
    src/bda/simplex_pan_recovery.c
    src/bda/simplex_pan_first_order.c
    src/bda/simplex_pan_cone_implicit.c
    src/bda/simplex_pan_basis.c
    src/bda/simplex_pan_basis_qr.c
    src/bda/simplex_pan_basis_spqr.c
    src/bda/simplex_pan_phase1.c
    src/bda/simplex_pan_standard.c
    src/bda/linalg/pan_dynamic_sparse_qr.c
    src/basis/simplex_basis.c
    src/basis/simplex_sparse_lu.c
    src/matrix/simplex_csc.c
    src/matrix/simplex_sparse_vector.c
    src/tableau/simplex_tableau_solver.c
    src/tableau/simplex_transform.c
    src/tableau/simplex_tableau.c
    src/tableau/simplex_phase.c
    src/tableau/simplex_pivot.c
    src/common/linalg.c
    src/common/utils.c)
add_library(lp_simplex SHARED ${SOURCES})
add_library(lp_simplex::lp_simplex ALIAS lp_simplex)
set_target_properties(lp_simplex PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
if(APPLE)
    set_target_properties(lp_simplex PROPERTIES MACOSX_RPATH ON INSTALL_NAME_DIR
                                                                "@rpath")
endif()
target_include_directories(
    lp_simplex PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(lp_simplex
                           PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})

target_link_libraries(lp_simplex PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(UNIX AND NOT APPLE)
    target_link_libraries(lp_simplex PRIVATE m)
endif()
target_compile_definitions(lp_simplex PRIVATE USE_BLAS USE_LAPACK)

set(LP_SIMPLEX_HAVE_GMP OFF)
set(LP_SIMPLEX_HAVE_SPQR OFF)
set(LP_SIMPLEX_HAVE_KLU OFF)

if(LP_SIMPLEX_USE_GMP)
    find_path(GMP_INCLUDE_DIR gmp.h)
    find_library(GMP_LIBRARY NAMES gmp)
    if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
        target_include_directories(lp_simplex SYSTEM PRIVATE ${GMP_INCLUDE_DIR})
        target_link_libraries(lp_simplex PRIVATE ${GMP_LIBRARY})
        set(LP_SIMPLEX_HAVE_GMP ON)
        target_compile_definitions(lp_simplex PRIVATE LP_SIMPLEX_HAVE_GMP)
        message(STATUS "Exact infeasibility verification: GMP")
    else()
        message(STATUS "GMP unavailable; exact basis verification disabled")
    endif()
endif()

if(LP_SIMPLEX_USE_SPQR)
    find_path(SPQR_INCLUDE_DIR SuiteSparseQR_C.h PATH_SUFFIXES suitesparse)
    find_library(SPQR_LIBRARY NAMES spqr)
    find_library(CHOLMOD_LIBRARY NAMES cholmod)
    find_library(SUITESPARSE_CONFIG_LIBRARY NAMES suitesparseconfig)
    if(SPQR_INCLUDE_DIR
       AND SPQR_LIBRARY
       AND CHOLMOD_LIBRARY
       AND SUITESPARSE_CONFIG_LIBRARY)
        target_include_directories(lp_simplex SYSTEM
                                   PRIVATE ${SPQR_INCLUDE_DIR})
        target_link_libraries(
            lp_simplex PRIVATE ${SPQR_LIBRARY} ${CHOLMOD_LIBRARY}
                               ${SUITESPARSE_CONFIG_LIBRARY})
        set(LP_SIMPLEX_HAVE_SPQR ON)
        target_compile_definitions(lp_simplex PRIVATE LP_SIMPLEX_HAVE_SPQR)
        target_sources(lp_simplex PRIVATE src/bda/linalg/pan_multifrontal_qr.c)
        message(STATUS "Using CHOLMOD/SPQR sparse Pan/BDA basis backend")
    else()
        message(
            STATUS "SuiteSparseQR not found; using Pan dense fallback backend")
    endif()
endif()

if(LP_SIMPLEX_USE_KLU)
    find_path(KLU_INCLUDE_DIR klu.h PATH_SUFFIXES suitesparse)
    find_library(KLU_LIBRARY NAMES klu)
    if(KLU_INCLUDE_DIR AND KLU_LIBRARY)
        target_include_directories(lp_simplex SYSTEM PRIVATE ${KLU_INCLUDE_DIR})
        target_link_libraries(lp_simplex PRIVATE ${KLU_LIBRARY})
        set(LP_SIMPLEX_HAVE_KLU ON)
        target_compile_definitions(lp_simplex PRIVATE LP_SIMPLEX_HAVE_KLU)
        message(STATUS "Using SuiteSparse KLU basis backend")
    else()
        message(
            STATUS "SuiteSparse KLU not found; using internal sparse LU backend"
        )
    endif()
endif()

###############################################################################
# Installation of the Library
###############################################################################

install(
    TARGETS lp_simplex
    EXPORT lp_simplexTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/lp_simplex
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h")
set(LP_SIMPLEX_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/lp_simplex")
configure_package_config_file(
    cmake/lp_simplexConfig.cmake.in
    "${PROJECT_BINARY_DIR}/lp_simplexConfig.cmake"
    INSTALL_DESTINATION ${LP_SIMPLEX_INSTALL_CMAKEDIR})
install(
    EXPORT lp_simplexTargets
    NAMESPACE lp_simplex::
    DESTINATION ${LP_SIMPLEX_INSTALL_CMAKEDIR})
install(FILES "${PROJECT_BINARY_DIR}/lp_simplexConfig.cmake"
        DESTINATION ${LP_SIMPLEX_INSTALL_CMAKEDIR})
install(FILES LICENSE LICENSE.GPL-3.0
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/licenses/lp_simplex)

###############################################################################
# Build Mode
#
# Choose the build mode either "Release" or "Debug" by
# cmake -DCMAKE_BUILD_TYPE=Release/Debug ..
###############################################################################

## (Default) The "Release" mode is for installation.
if(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    # Keep CMake's toolchain-specific Release optimization level (normally -O3).
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall")
    ##
    ## The "Debug" mode is for unit testings.
elseif(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
endif()

if(BUILD_TESTING)
    add_subdirectory(test)
endif()

if(LP_SIMPLEX_BUILD_BENCHMARKS)
    add_subdirectory(data)
endif()

if(LP_SIMPLEX_BUILD_EXAMPLES)
    add_subdirectory(data/simple_examples)
endif()
