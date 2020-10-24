
PROGRAM main

USE mod_step

IMPLICIT NONE

INTEGER  (KIND = 4)                                                         ::&
  i1, j1, n                          
REAL     (KIND = 8)                                                         ::&
  PI, st, en                  
REAL     (KIND = 8), DIMENSION(:,:), ALLOCATABLE                            ::&
  r, d 

PI = 4.0D0*ATAN(1.0D0)
  
n = 2000
ALLOCATE(r(n,n))
ALLOCATE(d(n,n))

! A very inefficient loop
DO i1 = 1,n
  DO j1 = 1,n
    d(i1,j1) = ABS(SIN(i1*PI/n)*COS(j1*PI/n))
  END DO
END DO

st = OMP_GET_WTIME()

! CALL step_v0(n,d,r)
! CALL step_omp(n,d,r) 
! CALL step_v1(n,d,r) 
CALL step_v2(n,d,r) 

en = OMP_GET_WTIME()

WRITE(*,'(A,F6.2,A)') 'Time = ', en-st, 'seconds'

! ! When testing, run with n = 5 and uncomment below
! WRITE(*,*) 'Max error = ', MAXVAL(ABS(r(1,:) - &
!   (/ 0.58778525229247325D0, 0.47552825814757682D0, 0.47552825814757671D0,&
!   0.58778525229247325D0, 0.58778525229247325D0 /)))
! WRITE(*,*) 'Max error = ', MAXVAL(ABS(r(2,:) - &
!   (/ 0.95105651629515364D0, 0.58778525229247314D0, 0.58778525229247292D0,&
!   0.95105651629515364D0, 0.95105651629515364D0 /)))
! WRITE(*,*) 'Max error = ', MAXVAL(ABS(r(3,:) - &
!   (/ 0.95105651629515375D0, 0.58778525229247314D0, 0.58778525229247303D0,&
!   0.95105651629515375D0, 0.95105651629515375D0 /)))
! WRITE(*,*) 'Max error = ', MAXVAL(ABS(r(4,:) - &
!   (/ 0.58778525229247336D0, 0.47552825814757682D0, 0.47552825814757671D0,&
!   0.58778525229247336D0, 0.58778525229247336D0 /)))
! WRITE(*,*) 'Max error = ', MAXVAL(ABS(r(5,:)))

DEALLOCATE(r,d)

END PROGRAM main