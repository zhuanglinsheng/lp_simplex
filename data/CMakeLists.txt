add_executable(lp_solve_mps lp_solve_mps.c)
target_link_libraries(lp_solve_mps PRIVATE lp_simplex)
if(UNIX AND NOT APPLE)
    target_link_libraries(lp_solve_mps PRIVATE m)
endif()

if(BUILD_TESTING)
    add_test(
        NAME test_mps_driver
        COMMAND
            ${CMAKE_COMMAND} -DDRIVER=$<TARGET_FILE:lp_solve_mps>
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/mps-driver-test -P
            ${PROJECT_SOURCE_DIR}/test/test_mps_driver.cmake)
endif()

add_subdirectory(netlib)
