cmake_minimum_required(VERSION 3.12)
project(installed_lp_consumer LANGUAGES C)
include(CTest)
find_package(lp_simplex CONFIG REQUIRED PATHS "${PACKAGE_ROOT}" NO_DEFAULT_PATH)
if(NOT TARGET lp_simplex::lp_simplex)
    message(
        FATAL_ERROR "Installed package has no lp_simplex::lp_simplex target")
endif()
foreach(feature KLU SPQR GMP)
    if(NOT DEFINED lp_simplex_WITH_${feature})
        message(FATAL_ERROR "Missing resolved capability ${feature}")
    endif()
endforeach()
# Debug is the repository's ASAN testing configuration, not the release artifact.
if(PACKAGE_ASAN AND NOT MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
endif()
add_executable(consumer_c main.c)
set_target_properties(consumer_c PROPERTIES C_STANDARD 90 C_STANDARD_REQUIRED
                                                          ON)
target_link_libraries(consumer_c PRIVATE lp_simplex::lp_simplex)
add_test(NAME installed_c COMMAND consumer_c)
if(TEST_CXX)
    enable_language(CXX)
    add_executable(consumer_cpp main.cpp)
    target_link_libraries(consumer_cpp PRIVATE lp_simplex::lp_simplex)
    add_test(NAME installed_cpp COMMAND consumer_cpp)
endif()
file(
    GENERATE
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/imported-$<CONFIG>.txt"
    CONTENT "$<TARGET_FILE:lp_simplex::lp_simplex>\n")
