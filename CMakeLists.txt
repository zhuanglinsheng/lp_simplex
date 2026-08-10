
###############################################################################
# Setting
###############################################################################

cmake_minimum_required(VERSION 3.12)
project(lp_simplex LANGUAGES C)
include(CTest)
option(LP_SIMPLEX_BUILD_TOOLS "Build command-line diagnostic tools" ON)
option(LP_SIMPLEX_USE_KLU "Use SuiteSparse KLU for sparse basis factorizations" ON)
option(LP_SIMPLEX_USE_SPQR "Use CHOLMOD/SPQR sparse algebra for Pan/BDA" ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_EXTENSIONS OFF)
if(MSVC)
	# Force to always compile with W4
	if(CMAKE_C_FLAGS_DEBUG MATCHES "/W[0-4]")
		string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
	else()
		set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /W4")
	endif()
else()
	# Update if necessary
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wno-long-long -pedantic -fsanitize=address")
endif()
##
## Installation setting
##
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")


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
	${PROJECT_SOURCE_DIR}/src/presolve
	${PROJECT_SOURCE_DIR}/src/dual
	${PROJECT_SOURCE_DIR}/src/degeneracy/pan
	${PROJECT_SOURCE_DIR}/src/degeneracy/pan/linalg
	${PROJECT_SOURCE_DIR}/src/basis
	${PROJECT_SOURCE_DIR}/src/matrix
	${PROJECT_SOURCE_DIR}/src/tableau
	${PROJECT_SOURCE_DIR}/src/common
)

set(SOURCES
	src/core/model.c
	src/core/mps.c
	src/core/simplex.c
	src/core/simplex_problem.c
	src/presolve/simplex_presolve.c
	src/presolve/simplex_presolve_activity.c
	src/presolve/simplex_presolve_postsolve.c
	src/presolve/simplex_presolve_profile.c
	src/presolve/simplex_presolve_queue.c
	src/presolve/simplex_presolve_record.c
	src/presolve/simplex_presolve_rules.c
	src/presolve/simplex_presolve_substitution.c
	src/dual/simplex_dual.c
	src/dual/simplex_dual_bounds.c
	src/dual/simplex_dual_state.c
	src/dual/simplex_dual_feasibility.c
	src/dual/simplex_dual_pricing.c
	src/dual/simplex_singleton_dual.c
	src/dual/simplex_degeneracy.c
	src/degeneracy/pan/simplex_pan.c
	src/degeneracy/pan/simplex_pan_basis.c
	src/degeneracy/pan/simplex_pan_basis_spqr.c
	src/degeneracy/pan/simplex_pan_phase1.c
	src/degeneracy/pan/simplex_pan_standard.c
	src/degeneracy/pan/linalg/pan_dynamic_sparse_qr.c
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
	src/common/utils.c
)
add_library(lp_simplex SHARED
	${SOURCES}
)
target_include_directories(lp_simplex PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)
target_include_directories(lp_simplex PRIVATE ${LP_SIMPLEX_PRIVATE_INCLUDE_DIRS})

target_link_libraries(lp_simplex PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_compile_definitions(lp_simplex PRIVATE USE_BLAS USE_LAPACK)

if(LP_SIMPLEX_USE_SPQR)
	find_path(SPQR_INCLUDE_DIR SuiteSparseQR_C.h PATH_SUFFIXES suitesparse)
	find_library(SPQR_LIBRARY NAMES spqr)
	find_library(CHOLMOD_LIBRARY NAMES cholmod)
	find_library(SUITESPARSE_CONFIG_LIBRARY NAMES suitesparseconfig)
	if(SPQR_INCLUDE_DIR AND SPQR_LIBRARY AND CHOLMOD_LIBRARY AND
	   SUITESPARSE_CONFIG_LIBRARY)
		target_include_directories(lp_simplex SYSTEM PRIVATE ${SPQR_INCLUDE_DIR})
		target_link_libraries(lp_simplex PRIVATE ${SPQR_LIBRARY}
			${CHOLMOD_LIBRARY} ${SUITESPARSE_CONFIG_LIBRARY})
		target_compile_definitions(lp_simplex PRIVATE LP_SIMPLEX_HAVE_SPQR)
		target_sources(lp_simplex PRIVATE
			src/degeneracy/pan/linalg/pan_multifrontal_qr.c)
		message(STATUS "Using CHOLMOD/SPQR sparse Pan/BDA basis backend")
	else()
		message(STATUS "SuiteSparseQR not found; using Pan dense fallback backend")
	endif()
endif()

if(LP_SIMPLEX_USE_KLU)
	find_path(KLU_INCLUDE_DIR klu.h PATH_SUFFIXES suitesparse)
	find_library(KLU_LIBRARY NAMES klu)
	if(KLU_INCLUDE_DIR AND KLU_LIBRARY)
		target_include_directories(lp_simplex SYSTEM PRIVATE ${KLU_INCLUDE_DIR})
		target_link_libraries(lp_simplex PRIVATE ${KLU_LIBRARY})
		target_compile_definitions(lp_simplex PRIVATE LP_SIMPLEX_HAVE_KLU)
		message(STATUS "Using SuiteSparse KLU basis backend")
	else()
		message(STATUS "SuiteSparse KLU not found; using internal sparse LU backend")
	endif()
endif()


###############################################################################
# Installation of the Library
###############################################################################

install(TARGETS lp_simplex LIBRARY
	DESTINATION lib
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/lp_simplex
	DESTINATION include
)


###############################################################################
# Build Mode
#
# Choose the build mode either "Release" or "Debug" by
# cmake -DCMAKE_BUILD_TYPE=Release/Debug ..
###############################################################################

## (Default) The "Release" mode is for installation.
if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
	# Keep CMake's toolchain-specific Release optimization level (normally -O3).
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall")
##
## The "Debug" mode is for unit testings.
elseif (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
endif ()

if(BUILD_TESTING)
	add_subdirectory(test)
endif()

if(LP_SIMPLEX_BUILD_TOOLS)
	add_subdirectory(tools)
endif()
