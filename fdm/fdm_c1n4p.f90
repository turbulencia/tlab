#include "types.h"

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C1N4P_LHS(imax, dx, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i

! ###################################################################
  DO i = 1,imax
     a(i) = C_025_R
     b(i) = C_1_R
     c(i) = C_025_R
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
  c(imax) = c(imax)*dx(1)
  b(1)    = b(1)   *dx(1)
  a(2)    = a(2)   *dx(1)

  DO i = 2,imax-1
     c(i-1) = c(i-1)*dx(i)
     b(i)   = b(i)  *dx(i)
     a(i+1) = a(i+1)*dx(i)
  ENDDO

  c(imax-1) = c(imax-1)*dx(imax)
  b(imax)   = b(imax)  *dx(imax)
  a(1)      = a(1)     *dx(imax)

  RETURN
END SUBROUTINE FDM_C1N4P_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N4P_RHS(imax,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk, srt,end,siz

! ###################################################################
!$omp parallel default(none) &
!$omp private(srt,end,siz,i,jk) &
!$omp shared(d,u,imax,jkmax) 
  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz) 

  DO jk = srt,end
     d(jk,1)    = C_075_R * (u(jk,2) - u(jk,imax))
     d(jk,imax) = C_075_R * (u(jk,1) - u(jk,imax-1))
  ENDDO

  DO i = 2,imax-1
     DO jk = srt,end
        d(jk,i) = C_075_R * (u(jk,i+1) - u(jk,i-1))
     ENDDO
  ENDDO

!$omp end parallel

  RETURN
END SUBROUTINE FDM_C1N4P_RHS
