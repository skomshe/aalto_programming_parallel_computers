
MODULE mod_step

    USE OMP_LIB

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE step_v0(n,d,r)

        IMPLICIT NONE

        INTEGER  (KIND = 4)                  :: n                         
        REAL     (KIND = 8) , DIMENSION(n,n) :: r, d 

        INTEGER  (KIND = 4)                  :: i1, j1

        DO i1 = 1,n
            DO j1 = 1,n
                r(i1,j1) = MINVAL(d(i1,1:n) + d(1:n,j1))
            END DO
        END DO
    
    END SUBROUTINE step_v0

    SUBROUTINE step_v1(n,d,r)

        IMPLICIT NONE

        INTEGER  (KIND = 4)                  :: n                         
        REAL     (KIND = 8) , DIMENSION(n,n) :: r, d 

        INTEGER  (KIND = 4)                  :: i1, j1
        REAL     (KIND = 8)                  :: tmp

        DO i1 = 1,n
            DO j1 = 1,n
                tmp = MINVAL(d(i1,1:n) + d(1:n,j1))
                r(i1,j1) = tmp
            END DO
        END DO
    
    END SUBROUTINE step_v1


END MODULE mod_step