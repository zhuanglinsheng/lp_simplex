NAME          VARBOUND
* Upper-only and multiple two-sided variable bounds
* EXPECTED STATUS: OPTIMAL
* EXPECTED OBJECTIVE: -5
ROWS
 N  OBJECTIV
 G  C1
COLUMNS
    X1        OBJECTIV             1
    X1        C1                   1
    X2        OBJECTIV             1
    X3        OBJECTIV            -1
RHS
    RHS1      C1                  -4
BOUNDS
 MI BND1      X1
 UP BND1      X1                   2
 LO BND1      X2                   3
 UP BND1      X2                   5
 LO BND1      X3                  -2
 UP BND1      X3                   4
ENDATA
