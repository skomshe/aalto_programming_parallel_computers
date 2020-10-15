
PROGRAM main

    USE mod_step
    
    IMPLICIT NONE
    
    INTEGER  (KIND = 4)                                 :: i1, n                          
    REAL     (KIND = 8)                                 :: st, en
    REAL     (KIND = 8) , DIMENSION(:,:), ALLOCATABLE   :: r, d 

    n = 2000 
    ALLOCATE(r(n,n))
    ALLOCATE(d(n,n))

    d = 1.0D0
    DO i1 = 1,n
        d(i1,i1) = 0.0D0
    END DO

    st = OMP_GET_WTIME()

    ! CALL step_v0(n,d,r)
    CALL step_v1(n,d,r) 

    en = OMP_GET_WTIME()

    WRITE(*,'(A,E15.4,A)') 'Time = ', en-st, 'seconds'

    DEALLOCATE(r,d)

END PROGRAM main