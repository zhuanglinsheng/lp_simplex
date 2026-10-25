foreach(required BUILD_DIR WORK_DIR CONSUMER_SOURCE CTEST_COMMAND CONFIG)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "Missing ${required}")
    endif()
endforeach()
# Keep cleanup confined to this CTest fixture, never a caller's install prefix.
if(NOT WORK_DIR STREQUAL "${BUILD_DIR}/test/installed-package")
    message(FATAL_ERROR "Unexpected package fixture directory: ${WORK_DIR}")
endif()
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
function(checked label)
    execute_process(
        COMMAND ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error)
    file(APPEND "${WORK_DIR}/steps.log" "${label}\n${output}\n${error}\n")
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "${label} failed (${result})\n${output}\n${error}")
    endif()
endfunction()
set(original "${WORK_DIR}/stage original")
set(relocated "${WORK_DIR}/stage relocated")
checked(
    install
    "${CMAKE_COMMAND}"
    -E
    env
    "DESTDIR="
    "${CMAKE_COMMAND}"
    "-DCMAKE_INSTALL_PREFIX=${original}"
    "-DCMAKE_INSTALL_CONFIG_NAME=${CONFIG}"
    -P
    "${BUILD_DIR}/cmake_install.cmake")
file(RENAME "${original}" "${relocated}")
if(EXISTS "${original}")
    message(FATAL_ERROR "Original prefix still exists after relocation")
endif()
foreach(license LICENSE LICENSE.GPL-3.0)
    if(NOT EXISTS "${relocated}/${DATAROOT}/licenses/lp_simplex/${license}")
        message(FATAL_ERROR "Installed package omitted ${license}")
    endif()
endforeach()
file(GLOB_RECURSE configs "${relocated}/*.cmake")
foreach(config ${configs})
    file(READ "${config}" content)
    foreach(forbidden "${original}" "${BUILD_DIR}" "${PROJECT_SOURCE}")
        string(FIND "${content}" "${forbidden}" found)
        if(NOT found EQUAL -1)
            message(
                FATAL_ERROR
                    "Export retains original/source/build path: ${config}")
        endif()
    endforeach()
endforeach()
set(extra)
if(DEFINED LIBDIR AND NOT LIBDIR STREQUAL "lib")
    list(APPEND extra
         "-Dlp_simplex_DIR=${relocated}/${LIBDIR}/cmake/lp_simplex")
endif()
if(GENERATOR_PLATFORM)
    list(APPEND extra -A "${GENERATOR_PLATFORM}")
endif()
if(GENERATOR_TOOLSET)
    list(APPEND extra -T "${GENERATOR_TOOLSET}")
endif()
if(CONFIG STREQUAL "Debug" AND DEBUG_ASAN)
    set(asan ON)
else()
    set(asan OFF)
endif()
checked(
    configure-consumer
    "${CMAKE_COMMAND}"
    -G
    "${GENERATOR}"
    ${extra}
    "-H${CONSUMER_SOURCE}"
    "-B${WORK_DIR}/consumer"
    "-DPACKAGE_ROOT=${relocated}"
    "-DCMAKE_BUILD_TYPE=${CONFIG}"
    "-DCMAKE_C_COMPILER=${C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
    "-DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF"
    "-DCMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY=ON"
    "-DTEST_CXX=${TEST_CXX}"
    "-DPACKAGE_ASAN=${asan}")
checked(build-consumer "${CMAKE_COMMAND}" --build "${WORK_DIR}/consumer"
        --config "${CONFIG}")
file(READ "${WORK_DIR}/consumer/imported-${CONFIG}.txt" imported)
string(FIND "${imported}" "${relocated}/" relocated_location)
if(NOT relocated_location EQUAL 0)
    message(
        FATAL_ERROR "Consumer did not import the relocated library: ${imported}"
    )
endif()
if(WIN32)
    set(ENV{PATH} "${relocated}/${BINDIR};$ENV{PATH}")
endif()
execute_process(
    COMMAND "${CTEST_COMMAND}" -C "${CONFIG}" --output-on-failure
    WORKING_DIRECTORY "${WORK_DIR}/consumer"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error)
file(APPEND "${WORK_DIR}/steps.log" "test-consumer\n${output}\n${error}\n")
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Installed consumer failed\n${output}\n${error}")
endif()
message(STATUS "Relocated installed package and C/C++ consumers passed")
