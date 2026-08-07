
###############################################################################
# Setting
###############################################################################

cmake_minimum_required(VERSION 3.12)
project(lp_simplex LANGUAGES C)
include(CTest)
option(LP_SIMPLEX_BUILD_TOOLS "Build command-line diagnostic tools" ON)
option(LP_SIMPLEX_USE_KLU "Use SuiteSparse KLU for sparse basis factorizations" ON)
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

set(SOURCES
	src/model.c
	src/mps.c
	src/simplex.c
	src/simplex_tableau_solver.c
	src/simplex_dual.c
	src/simplex_singleton_dual.c
	src/simplex_csc.c
	src/simplex_basis.c
	src/simplex_sparse_lu.c
	src/simplex_transform.c
	src/simplex_tableau.c
	src/simplex_phase.c
	src/simplex_pivot.c
	src/linalg.c
	src/utils.c
)
add_library(lp_simplex SHARED
	${SOURCES}
)
target_include_directories(lp_simplex PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)
target_include_directories(lp_simplex PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(lp_simplex PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_compile_definitions(lp_simplex PRIVATE USE_BLAS USE_LAPACK)

if(LP_SIMPLEX_USE_KLU)
	find_path(KLU_INCLUDE_DIR klu.h PATH_SUFFIXES suitesparse)
	find_library(KLU_LIBRARY NAMES klu)
	if(KLU_INCLUDE_DIR AND KLU_LIBRARY)
		target_include_directories(lp_simplex PRIVATE ${KLU_INCLUDE_DIR})
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
