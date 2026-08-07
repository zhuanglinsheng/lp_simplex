NAME          FREELOW
* Free variable combined with a negative lower bound
* EXPECTED STATUS: OPTIMAL
* EXPECTED OBJECTIVE: -22
ROWS
 N  OBJECTIV
 L  C1
 G  C2
COLUMNS
    X1        OBJECTIV            -1
    X1        C1                  -3
    X1        C2                  -1
    X2        OBJECTIV             4
    X2        C1                   1
    X2        C2                  -2
RHS
    RHS1      C1                   6
    RHS1      C2                  -4
BOUNDS
 FR BND1      X1
 LO BND1      X2                  -3
ENDATA
