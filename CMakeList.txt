###############################################################################
# Setting
###############################################################################

cmake_minimum_required(VERSION 3.12)
project(lp_simplex LANGUAGES C)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_EXTENSIONS OFF)
if(MSVC)
	# Force to always compile with W4
	if(CMAKE_C_FLAGS_DEBUG MATCHES "/W[0-4]")
		string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
	else()
		set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /W4 -fsanitize=address")
	endif()
else()
	# Update if necessary
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wno-long-long -pedantic -fsanitize=address")
endif()
##
## Installation setting
##
## If user would like to choose their own default installation location,
## set the value of `CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT` to be
## `true` by switch the comments:
##
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "..." FORCE)
set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT false)
##
## Then set the installation location below by changing the default value
## of `CMAKE_INSTALL_PREFIX`:
##
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "$HOME/.local" CACHE PATH "..." FORCE)
endif()
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
##
## Third-party libraries
## add_compile_definitions(WITH_BLAS)


###############################################################################
# shared library
###############################################################################

set(SOURCES
	src/lp_simplex_bsc.c
	src/lp_simplex_std.c
	src/lp_simplex_gen.c
)
add_library(lp_simplex SHARED
	${SOURCES}
)
target_include_directories(lp_simplex PRIVATE
	${PROJECT_SOURCE_DIR}/include
)

find_package(BLAS REQUIRED)

if(BLAS_FOUND)
	message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
	set(USE_BLAS TRUE)
else()
	message(WARNING "BLAS not found")
	set(USE_BLAS FALSE)
endif()

if(USE_BLAS)
	target_link_libraries(lp_simplex ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
	target_compile_definitions(lp_simplex PRIVATE USE_BLAS USE_LAPACK)
else()
	target_link_libraries(lp_simplex)
endif()


###############################################################################
# Installation
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
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall -O2 -DNDEBUG")
##
## The "Debug" mode is for unit testings.
elseif (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
	enable_testing()
	add_subdirectory(test)
endif ()
