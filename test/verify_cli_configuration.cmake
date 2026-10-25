# Check both the process result and the actual solver/driver mode. A CTest
# PASS_REGULAR_EXPRESSION by itself ignores the process exit code.
execute_process(
    COMMAND "${DRIVER}" ${DRIVER_ARGUMENTS} --expected "${EXPECTED}" "${MODEL}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE diagnostics)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "Driver failed (${result}):\n${output}\n${diagnostics}")
endif()
if(NOT output MATCHES "configuration: experimental-environment")
    message(FATAL_ERROR "Missing experimental driver mode:\n${output}")
endif()
if(NOT diagnostics MATCHES "execution-config: mode=experimental-environment")
    message(
        FATAL_ERROR
            "Solver did not enter the experimental root:\n${diagnostics}")
endif()
