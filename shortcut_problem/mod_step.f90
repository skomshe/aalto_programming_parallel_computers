
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

SUBROUTINE step_v2(n,na,nb,nab,d,r)

IMPLICIT NONE

INTEGER  (KIND = 4)                                                         ::&
  n, na, nb, nab
REAL     (KIND = 8), DIMENSION(:,:)                                         ::&
  r, d 

INTEGER  (KIND = 4)                                                         ::&
  i1, j1, k1, idx
INTEGER  (KIND = 4), DIMENSION(nb)                                          ::&
  idxA
REAL     (KIND = 8)                                                         ::&
  tmp
REAL     (KIND = 8), DIMENSION(nb)                                          ::&
  tmpA, vv
REAL     (KIND = 8), DIMENSION(nab,n)                                       ::&
  d_pad
REAL     (KIND = 8), DIMENSION(n,nab)                                       ::&
  t_pad

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i1,j1,tmp) &
!$OMP SHARED(n,d,d_pad,t_pad)   
DO i1 = 1,n
  DO j1 = 1,n
    tmp = d(j1,i1)
    d_pad(j1,i1) = tmp
    t_pad(i1,j1) = tmp
  END DO
END DO
!$OMP END PARALLEL DO
d_pad(n+1:nab,1:n) = 1.0D5
t_pad(1:n,n+1:nab) = 1.0D5

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i1,j1,k1,idx,idxA,tmpA,vv) &
!$OMP SHARED(n,na,nb,r,d_pad,t_pad)        
DO i1 = 1,n
  DO j1 = 1,n

    vv(1) = t_pad(1,i1) + d_pad(1,j1)
    vv(2) = t_pad(2,i1) + d_pad(2,j1)
    vv(3) = t_pad(3,i1) + d_pad(3,j1)
    vv(4) = t_pad(4,i1) + d_pad(4,j1)

    DO k1 = 2,na
      idx = (k1-1)*nb

      idxA(1) = idx + 1
      idxA(2) = idx + 2
      idxA(3) = idx + 3
      idxA(4) = idx + 4

      tmpA(1) = t_pad(idxA(1),i1) + d_pad(idxA(1),j1)
      tmpA(2) = t_pad(idxA(2),i1) + d_pad(idxA(2),j1)
      tmpA(3) = t_pad(idxA(3),i1) + d_pad(idxA(3),j1)
      tmpA(4) = t_pad(idxA(4),i1) + d_pad(idxA(4),j1)

      vv(1) = MIN(tmpA(1),vv(1))
      vv(2) = MIN(tmpA(2),vv(2))
      vv(3) = MIN(tmpA(3),vv(3))
      vv(4) = MIN(tmpA(4),vv(4))
    END DO

    r(i1,j1) = MINVAL(vv)
    
  END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE step_v2
!------------------------------------------------------------------------------


END MODULE mod_step