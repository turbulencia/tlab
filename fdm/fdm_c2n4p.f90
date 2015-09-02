#include "types.h"

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C2N4P_LHS(imax, dx, a,b,c)

  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i

! #######################################################################
  DO i = 1,imax
     a(i) = C_01_R
     b(i) = C_1_R
     c(i) = C_01_R
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
  c(imax) = c(imax)*(dx(1)**2)
  b(1)    = b(1)   *(dx(1)**2)
  a(2)    = a(2)   *(dx(1)**2)

  DO i = 2,imax-1
     c(i-1) = c(i-1)*(dx(i)**2)
     b(i)   = b(i)  *(dx(i)**2)
     a(i+1) = a(i+1)*(dx(i)**2)
  ENDDO

  c(imax-1) = c(imax-1)*(dx(imax)**2)
  b(imax)   = b(imax)  *(dx(imax)**2)
  a(1)      = a(1)     *(dx(imax)**2)

  RETURN
END SUBROUTINE FDM_C2N4P_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C2N4P_RHS(imax,jkmax, u,d)

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk
  TREAL c65dx2

! #######################################################################
!$omp parallel default(none) &
!$omp private(c65dx2,jk,i) & 
!$omp shared(u,d,jkmax,imax) 

  c65dx2 = C_6_R/C_5_R

  DO jk = 1,jkmax
     d(jk,1)    = c65dx2*(u(jk,2)-C_2_R*u(jk,1)   +u(jk,imax)  )
     d(jk,imax) = c65dx2*(u(jk,1)-C_2_R*u(jk,imax)+u(jk,imax-1))
  ENDDO

  DO i = 2,imax-1
     DO jk = 1,jkmax
        d(jk,i) = c65dx2*(u(jk,i+1)-C_2_R*u(jk,i)+u(jk,i-1))
     ENDDO
  ENDDO
!$omp end parallel 
  RETURN
END SUBROUTINE FDM_C2N4P_RHS
