
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

SUBROUTINE step_v2(n,d,r)

IMPLICIT NONE

INTEGER  (KIND = 4)                                                         ::&
  n                         
REAL     (KIND = 8), DIMENSION(:,:)                                         ::&
  r, d 

INTEGER  (KIND = 4)                                                         ::&
  nb, na, nab, i1, j1, k1, k2, idx, idx1, idx2, idx3, idx4
REAL     (KIND = 8)                                                         ::&
  tmp, tmp1, tmp2, tmp3, tmp4
REAL     (KIND = 8), DIMENSION(:), ALLOCATABLE                              ::&
  vv
REAL     (KIND = 8), DIMENSION(:,:), ALLOCATABLE                            ::&
  d_pad, t_pad

nb = 4
na = (n + nb - 1)/nb
nab = na*nb

ALLOCATE(vv(nb))

ALLOCATE(d_pad(nab,n))
ALLOCATE(t_pad(n,nab))
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
!$OMP PRIVATE(i1,j1,k1,k2,idx,idx1,idx2,idx3,idx4,tmp,tmp1,tmp2,tmp3,tmp4,vv) &
!$OMP SHARED(n,na,nb,r,d_pad,t_pad)        
DO i1 = 1,n
  DO j1 = 1,n

    ! DO k1 = 1,nb
    !   vv(k1) = t_pad(k1,i1) + d_pad(k1,j1)
    ! END DO

    vv(1) = t_pad(1,i1) + d_pad(1,j1)
    vv(2) = t_pad(2,i1) + d_pad(2,j1)
    vv(3) = t_pad(3,i1) + d_pad(3,j1)
    vv(4) = t_pad(4,i1) + d_pad(4,j1)

    DO k1 = 2,na
      ! DO k2 = 1,nb
      !   idx = (k1-1)*nb + k2
      !   tmp = t_pad(idx,i1) + d_pad(idx,j1)
      !   vv(k2) = MIN(tmp,vv(k2))
      ! END DO

      idx = (k1-1)*nb

      idx1 = idx + 1
      idx2 = idx + 2
      idx3 = idx + 3
      idx4 = idx + 4

      tmp1 = t_pad(idx1,i1) + d_pad(idx1,j1)
      tmp2 = t_pad(idx2,i1) + d_pad(idx2,j1)
      tmp3 = t_pad(idx3,i1) + d_pad(idx3,j1)
      tmp4 = t_pad(idx4,i1) + d_pad(idx4,j1)

      vv(1) = MIN(tmp1,vv(1))
      vv(2) = MIN(tmp2,vv(2))
      vv(3) = MIN(tmp3,vv(3))
      vv(4) = MIN(tmp4,vv(4))

    END DO

    r(i1,j1) = MINVAL(vv)
    
  END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(vv,d_pad,t_pad)

END SUBROUTINE step_v2
!------------------------------------------------------------------------------


END MODULE mod_step