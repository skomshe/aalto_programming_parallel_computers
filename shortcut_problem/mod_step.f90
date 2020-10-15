
MODULE mod_step

USE omp_lib

IMPLICIT NONE

CONTAINS
!------------------------------------------------------------------------------

SUBROUTINE step_v0(n,d,r)

IMPLICIT NONE

INTEGER  (KIND = 4)                                                         ::&
  n                         
REAL     (KIND = 8), DIMENSION(:,:)                                         ::&
  r, d 

INTEGER  (KIND = 4)                                                         ::&
  i1, j1, k1
REAL     (KIND = 8)                                                         ::&
  v, z

DO i1 = 1,n
  DO j1 = 1,n
    v = d(i1,1) + d(1,j1)
    DO k1 = 2,n
      z = d(i1,k1) + d(k1,j1)
      v = MIN(v,z)
    END DO
    r(i1,j1) = v
  END DO
END DO

END SUBROUTINE step_v0
!------------------------------------------------------------------------------

SUBROUTINE step_omp(n,d,r)

IMPLICIT NONE

INTEGER  (KIND = 4)                                                         ::&
  n                         
REAL     (KIND = 8), DIMENSION(:,:)                                         ::&
  r, d 

INTEGER  (KIND = 4)                                                         ::&
  i1, j1, k1
REAL     (KIND = 8)                                                         ::&
  v, z

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i1,j1,k1,v,z) &
!$OMP SHARED(n,r,d)        
DO i1 = 1,n
  DO j1 = 1,n
    v = d(i1,1) + d(1,j1)
    DO k1 = 2,n
      z = d(i1,k1) + d(k1,j1)
      v = MIN(v,z)
    END DO
    r(i1,j1) = v
  END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE step_omp
!------------------------------------------------------------------------------

SUBROUTINE step_v1(n,d,r)

IMPLICIT NONE

INTEGER  (KIND = 4)                                                         ::&
  n                         
REAL     (KIND = 8), DIMENSION(:,:)                                         ::&
  r, d 

INTEGER  (KIND = 4)                                                         ::&
  i1, j1, k1
REAL     (KIND = 8)                                                         ::&
  v, z
REAL     (KIND = 8), DIMENSION(:,:), ALLOCATABLE                            ::&
  t

ALLOCATE(t(n,n))
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i1,j1) &
!$OMP SHARED(n,d,t)   
DO i1 = 1,n
  DO j1 = 1,n
    t(i1,j1) = d(j1,i1)
  END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i1,j1,k1,v,z) &
!$OMP SHARED(n,r,d,t)        
DO i1 = 1,n
  DO j1 = 1,n
    v = t(1,i1) + d(1,j1)
    DO k1 = 2,n
      z = t(k1,i1) + d(k1,j1)
      v = MIN(v,z)
    END DO
    r(i1,j1) = v
  END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(t)

END SUBROUTINE step_v1
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------


END MODULE mod_step