#include "types.h"

!########################################################################
!Left-hand side; tridiagonal matrix of the linear system
!########################################################################
SUBROUTINE FDM_C2N4_LHS(imax, imin_set_zero,imax_set_zero, dx,a,b,c)

  IMPLICIT NONE

  TINTEGER,                  INTENT(IN) :: imax, imin_set_zero, imax_set_zero
  TREAL,   DIMENSION(imax,2),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),  INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i
  TREAL vmult_imin, vmult_imax

!########################################################################
  vmult_imin  = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin  = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

!########################################################################
! third-order biased
  a(1)    = C_0_R
  b(1)    = C_1_R
  c(1)    = C_11_R * vmult_imin

  a(imax) = C_11_R * vmult_imax
  b(imax) = C_1_R
  c(imax) = C_0_R

  DO i=2, imax-1
     a(i) = C_01_R
     b(i) = C_1_R
     c(i) = C_01_R
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
  c(imax) = c(imax)*(dx(1,1)**2)
  b(1)    = b(1)   *(dx(1,1)**2)
  a(2)    = a(2)   *(dx(1,1)**2)

  DO i = 2,imax-1
     c(i-1) = c(i-1)*(dx(i,1)**2)
     b(i)   = b(i)  *(dx(i,1)**2)
     a(i+1) = a(i+1)*(dx(i,1)**2)
  ENDDO

  c(imax-1) = c(imax-1)*(dx(imax,1)**2)
  b(imax)   = b(imax)  *(dx(imax,1)**2)
  a(1)      = a(1)     *(dx(imax,1)**2)

  RETURN
END SUBROUTINE FDM_C2N4_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C2N4_RHS(iunif, imax,jkmax, imin_set_zero,imax_set_zero, dx, u,up,d)

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: iunif, imax,jkmax, imin_set_zero,imax_set_zero
  TREAL,   DIMENSION(imax,2),    INTENT(IN) :: dx
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u, up
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk
  TREAL vmult_imin, vmult_imax
  TREAL c65dx2, c1dx2, c13dx2, c27dx2, c15dx2

!########################################################################
  vmult_imin  = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin  = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

!########################################################################
! Left (i=1) boundary, third-order
  c1dx2  = C_1_R * vmult_imin
  c13dx2 = C_13_R * c1dx2
  c27dx2 = C_27_R * c1dx2
  c15dx2 = C_15_R * c1dx2

  DO jk = 1,jkmax
     d(jk,1) = c13dx2*u(jk,1) - c27dx2*u(jk,2) + c15dx2*u(jk,3) - c1dx2*u(jk,4)
  ENDDO

! Right (i=imax) boundary, third-order
  c1dx2  = C_1_R * vmult_imax
  c13dx2 = C_13_R * c1dx2
  c27dx2 = C_27_R * c1dx2
  c15dx2 = C_15_R * c1dx2

  DO jk = 1,jkmax
     d(jk,imax) = c13dx2*u(jk,imax) - c27dx2*u(jk,imax-1) + c15dx2*u(jk,imax-2) - c1dx2*u(jk,imax-3)
  ENDDO

! -------------------------------------------------------------------
! Uniform case
! -------------------------------------------------------------------
  IF ( iunif .EQ. 0 ) THEN
     c65dx2 = C_6_R/C_5_R
     DO i = 2,imax-1
        DO jk = 1,jkmax
           d(jk,i) = c65dx2*(u(jk,i+1)-C_2_R*u(jk,i)+u(jk,i-1))
        ENDDO
     ENDDO

! -------------------------------------------------------------------
! Nonuniform case
! -------------------------------------------------------------------
  ELSE
! Not yet developed
  ENDIF

  RETURN
END SUBROUTINE FDM_C2N4_RHS
