!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/10/01 - J.P. Mellado
!#              Fixed
!# 2010/03/09 - J.P. Mellado
!#              Splitting into two routines
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 8th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic.
!# Derived from FDM_C1N3D6.
!# First point with biased third-order    Eq. 4.1.3 with \alpha=2.
!# Second point with fourth-order scheme  Eq. 2.1.6 with \alpha=1/4.
!# Third point with sixth-order scheme    Eq. 2.1.7 (\alpha=1/3).
!# Interior point with eight-order scheme Eq. 2.1.13 with \alpha=3/8.
!#
!########################################################################
!# ARGUMENTS 
!#
!# u    In    function to be diferentiated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################
#include "types.h"

!########################################################################
!Left-hand side; tridiagonal matrix of the linear system
!########################################################################
SUBROUTINE FDM_C1N8_LHS(imax, imin_set_zero,imax_set_zero, dx,a,b,c)

  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax, imin_set_zero,imax_set_zero
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i
  TREAL vmult_imin, vmult_imax

! ###################################################################
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! ###################################################################
! Tridiagonal matrix of the linear system
! ###################################################################
! third-order biased
  a(1)    = C_0_R
  b(1)    = C_1_R
  c(1)    = C_2_R * vmult_imin

  a(imax) = C_2_R * vmult_imax
  b(imax) = C_1_R
  c(imax) = C_0_R

! fourth-order centered
  a(2)      = C_025_R
  b(2)      = C_1_R
  c(2)      = C_025_R

  a(imax-1) = C_025_R
  b(imax-1) = C_1_R
  c(imax-1) = C_025_R

! sixth-order centered
  a(3)      = C_1_R/C_3_R
  b(3)      = C_1_R
  c(3)      = C_1_R/C_3_R

  a(imax-2) = C_1_R/C_3_R
  b(imax-2) = C_1_R
  c(imax-2) = C_1_R/C_3_R

  DO i = 4,imax-3
     a(i) = C_3_R/C_8_R
     b(i) = C_1_R
     c(i) = C_3_R/C_8_R
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
  b(imax  ) = b(imax)  *dx(imax)
  a(1)      = a(1)     *dx(imax)

  RETURN
END SUBROUTINE FDM_C1N8_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N8_RHS(imax,jkmax, imin_set_zero,imax_set_zero, u,d)
#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax, imin_set_zero, imax_set_zero
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk
  TINTEGER srt,end,siz
  TREAL vmult_imin, vmult_imax, c18dx, c36dx, c34dx
  TREAL c52dx1, c2dx1, c12dx1, c52dxmx, c2dxmx, c12dxmx
  TREAL c32dx, c20dx, c48dx

! #######################################################################
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

!$omp parallel default( none ) &
!$omp private( jk,c12dx1,c52dx1,c2dx1,c12dxmx,c52dxmx,c2dxmx,c34dx,c18dx,c36dx,c32dx,c20dx,c48dx, srt,end,siz ) &
!$omp shared(d,u,imax,vmult_imin,vmult_imax,jkmax)
! #######################################################################

  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz) 

  c12dx1 = C_05_R * vmult_imin
  c52dx1 = C_5_R  * c12dx1
  c2dx1  = C_2_R  * vmult_imin

  c12dxmx = C_05_R * vmult_imax
  c52dxmx = C_5_R  * c12dxmx
  c2dxmx  = C_2_R  * vmult_imax

  c34dx = C_075_R
  c18dx = C_14_R/C_18_R
  c36dx = C_1_R /C_36_R
  DO jk=srt,end
! third-order
     d(jk,1)    =-c52dx1 *u(jk,1)    + c2dx1 *u(jk,2)      + c12dx1 *u(jk,3)
     d(jk,imax) = c52dxmx*u(jk,imax) - c2dxmx*u(jk,imax-1) - c12dxmx*u(jk,imax-2)
! fourth-order centered
     d(jk,2)      = c34dx*(u(jk,3)    - u(jk,1))
     d(jk,imax-1) = c34dx*(u(jk,imax) - u(jk,imax-2))
! sixth-order centered
     d(jk,3)      = c18dx*(u(jk,4)      - u(jk,2))      + c36dx*(u(jk,5)    - u(jk,1))
     d(jk,imax-2) = c18dx*(u(jk,imax-1) - u(jk,imax-3)) + c36dx*(u(jk,imax) - u(jk,imax-4))
  ENDDO

  c32dx = C_5_R*C_5_R/(C_4_R*C_8_R)
  c20dx = C_1_R      /(C_2_R*C_10_R)
  c48dx =-C_1_R      /(C_3_R*C_10_R*C_16_R)
  DO i = 4,imax-3
     DO jk = srt,end
        d(jk,i) = c32dx*(u(jk,i+1) - u(jk,i-1)) + c20dx*(u(jk,i+2) - u(jk,i-2)) &
                + c48dx*(u(jk,i+3) - u(jk,i-3))
     ENDDO
  ENDDO
!$omp end parallel

  RETURN
END SUBROUTINE FDM_C1N8_RHS

