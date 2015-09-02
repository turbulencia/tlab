#include "types.h"

!########################################################################
!Left-hand side; tridiagonal matrix of the linear system
!########################################################################
SUBROUTINE FDM_C1N4_LHS(imax, imin_set_zero,imax_set_zero, dx,a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax, imin_set_zero, imax_set_zero
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i
  TREAL vmult_imin, vmult_imax

! #######################################################################
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! #######################################################################
! third-order biased
  a(1)    = C_0_R
  b(1)    = C_1_R
  c(1)    = C_2_R * vmult_imin

  a(imax) = C_2_R * vmult_imax
  b(imax) = C_1_R
  c(imax) = C_0_R

  DO i = 2,imax-1
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
     b(i)   = b(i  )*dx(i)
     a(i+1) = a(i+1)*dx(i)
  ENDDO

  c(imax-1) = c(imax-1)*dx(imax)
  b(imax)   = b(imax)  *dx(imax)
  a(1)      = a(1)     *dx(imax)

  RETURN
END SUBROUTINE FDM_C1N4_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N4_RHS(imax,jkmax, imin_set_zero,imax_set_zero, u,d)
  
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) ::imax, jkmax, imin_set_zero, imax_set_zero
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk, srt,end,siz
  TREAL vmult_imin, vmult_imax, c34dx
  TREAL c52dx1, c2dx1, c12dx1, c52dxmx, c2dxmx, c12dxmx

! #######################################################################
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

!$omp parallel default(none) &
!$omp private (c12dx1, c52dx1,c2dx1,c12dxmx,c52dxmx,c2dxmx,c34dx,srt,end,siz) &
!$omp shared(u,d,jkmax,imax,vmult_imin,vmult_imax) 

  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz) 

! #######################################################################
! third-order biased
  c12dx1 = C_05_R * vmult_imin
  c52dx1 = C_5_R  * c12dx1
  c2dx1  = C_2_R  * vmult_imin

  c12dxmx = C_05_R * vmult_imax
  c52dxmx = C_5_R  * c12dxmx
  c2dxmx  = C_2_R  * vmult_imax

  DO jk = srt,end
     d(jk,1) = -c52dx1*u(jk,1) + c2dx1*u(jk,2) + c12dx1*u(jk,3)
     d(jk,imax) = c52dxmx*u(jk,imax) - c2dxmx*u(jk,imax-1) - c12dxmx*u(jk,imax-2)
  ENDDO

  c34dx = C_075_R
  DO i = 2,imax-1
     DO jk = srt,end
        d(jk,i) = c34dx * (u(jk,i+1) - u(jk,i-1))
     ENDDO
  ENDDO
!$omp end parallel 
  RETURN
END SUBROUTINE FDM_C1N4_RHS
