
PROGRAM main

    USE mod_step
    USE OMP_LIB
    
    IMPLICIT NONE
    
    INTEGER  (KIND = 4)                  :: i1
    INTEGER  (KIND = 4) , PARAMETER      :: n = 2000                           
    REAL     (KIND = 8)                  :: st, en
    REAL     (KIND = 8) , DIMENSION(n,n) :: r, d 

    CALL RANDOM_NUMBER(d)
    DO i1 = 1,n
        d(i1,i1) = 0.0D0
    END DO

    CALL CPU_TIME(st)

    ! CALL step_v0(n,d,r)
    ! CALL step_v1(n,d,r)

    CALL CPU_TIME(en)

    WRITE(*,'(A,E15.4,A)') 'Time = ', en-st, 'seconds'

END PROGRAM main