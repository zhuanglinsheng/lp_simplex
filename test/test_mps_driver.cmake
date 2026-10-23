# A generic driver must not infer expectations from a dataset name or path.
file(MAKE_DIRECTORY "${WORK_DIR}/infeasible")
set(model "${WORK_DIR}/infeasible/afiro.mps")
set(contents "NAME          AFIRO
ROWS
 N  COST
 G  R1
COLUMNS
    X         COST                 1
    X         R1                   1
RHS
    RHS1      R1                   1
ENDATA
")
file(WRITE "${model}" "${contents}")

function(check_result expected_code expected_text)
    execute_process(COMMAND "${DRIVER}" ${ARGN} "${model}"
        RESULT_VARIABLE code OUTPUT_VARIABLE output ERROR_VARIABLE error
        TIMEOUT 15)
    if(NOT "${code}" STREQUAL "${expected_code}" OR
       NOT output MATCHES "${expected_text}")
        message(FATAL_ERROR "Driver returned ${code}:\n${output}\n${error}")
    endif()
endfunction()

check_result(0 "no reference; successful solve required")
check_result(0 "prediction reached" --expected 1)
check_result(1 "prediction not reached" --expected 2)
check_result(1 "prediction not reached" --infeasible)

# A genuinely infeasible model passes only with an explicit expectation.
string(REPLACE "ENDATA" "BOUNDS\n UP BND       X                    0\nENDATA" contents "${contents}")
file(WRITE "${model}" "${contents}")
check_result(0 "prediction reached" --infeasible)
check_result(1 "prediction not reached")

# Progress is separate from the stable result output, and a feasible point is
# checkpointed even without verbose progress. Invalid intervals are rejected.
file(WRITE "${model}" "NAME          PROGRESS
ROWS
 N  COST
 G  R1
COLUMNS
    X         COST                 1
    X         R1                   1
RHS
    RHS1      R1                   1
ENDATA
")
set(point "${WORK_DIR}/incumbent.point")
execute_process(COMMAND "${DRIVER}" --progress --save-point "${point}" "${model}"
    RESULT_VARIABLE code OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 15)
if(NOT "${code}" STREQUAL "0" OR NOT error MATCHES "progress:.*lower=" OR
   NOT error MATCHES "upper=.*gap=" OR NOT error MATCHES "stage=finish" OR
   NOT EXISTS "${point}")
    message(FATAL_ERROR "Missing progress/checkpoint: ${output}\n${error}")
endif()
file(READ "${point}" coordinates)
if(NOT coordinates STREQUAL "1\n" OR EXISTS "${point}.tmp")
    message(FATAL_ERROR "Incomplete checkpoint: ${coordinates}")
endif()
check_result(2 "Usage:" --progress-interval 0)
check_result(2 "Usage:" --progress-interval nan)
check_result(2 "Usage:" --progress-interval)

# A failed new write cannot destroy the last complete checkpoint.
file(WRITE "${point}" "previous complete checkpoint\n")
file(MAKE_DIRECTORY "${point}.tmp")
execute_process(COMMAND "${DRIVER}" --save-point "${point}" "${model}"
    RESULT_VARIABLE code OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 15)
file(READ "${point}" coordinates)
if(NOT "${code}" STREQUAL "1" OR
   NOT coordinates STREQUAL "previous complete checkpoint\n" OR
   NOT error MATCHES "could not checkpoint")
    message(FATAL_ERROR "Failed checkpoint destroyed its predecessor")
endif()
file(REMOVE_RECURSE "${point}.tmp")
