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
!# 2009/12/09 - J.P. Mellado
!#              Splitting into two routines
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic.
!# Interior points according to Eq. 2.1.7 (\alpha=1/3).
!# Carpenter, Gottlieb and Aberbanel, JCP, 1993, study the effect of 
!# boundary points on stability. However, the scheme 3-4-6-4-3 implemented
!# here seems to work better than the 3-5-6-5-3 proposed there. This
!# fourth-order scheme is Eq. 2.1.6 with \alpha=1/4.
!# The BCs are imposed with biased third-order Eq. 4.1.3 with \alpha=2.
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
SUBROUTINE FDM_C1N6_LHS(imax, imin_set_zero,imax_set_zero, dx,a,b,c)

  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax, imin_set_zero,imax_set_zero
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i
  TREAL vmult_imin, vmult_imax, c0103

! #######################################################################
  c0103 = C_1_R/C_3_R

  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! #######################################################################
! third-order biased
  a(1)    = C_0_R
  b(1)    = C_1_R
  c(1)    = C_2_R *vmult_imin

  a(imax) = C_2_R *vmult_imax
  b(imax) = C_1_R
  c(imax) = C_0_R

! fourth-order centered
!   a(2)      = C_025_R
!   b(2)      = C_1_R
!   c(2)      = C_025_R

!   a(imax-1) = C_025_R
!   b(imax-1) = C_1_R
!   c(imax-1) = C_025_R
! fifth-order biased
  a(2) = C_1_R/C_6_R
  b(2) = C_1_R
  c(2) = C_1_R/C_2_R
  
  a(imax-1) = C_1_R/C_2_R
  b(imax-1) = C_1_R
  c(imax-1) = C_1_R/C_6_R

  DO i = 3,imax-2
     a(i) = c0103
     b(i) = C_1_R
     c(i) = c0103
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
END SUBROUTINE FDM_C1N6_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N6_RHS(imax,jkmax, imin_set_zero,imax_set_zero, u,d)
#ifdef USE_OPENMP
  USE OMP_LIB
#endif 

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax,jkmax, imin_set_zero,imax_set_zero
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d
  TINTEGER                                  :: srt,end,siz    ! OpenMP Partition variables

! -------------------------------------------------------------------
  TINTEGER i, jk
  TREAL vmult_imin, vmult_imax
  TREAL c1418, c0136, c0118, c0102, c0509 !, c0304
  TREAL c52dx1, c2dx1, c12dx1, c52dxmx, c2dxmx, c12dxmx

! #######################################################################
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! #######################################################################
!$omp parallel default ( none ) &
!$omp private( jk,i, c1418,c0136,c0118,c0102,c0509, c12dx1,c52dx1,c2dx1, c12dxmx,c52dxmx,c2dxmx, srt,end,siz ) &
!$omp shared(d,u,jkmax,vmult_imin,vmult_imax,imax) 

  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz)

  c12dx1  = C_05_R  *vmult_imin
  c52dx1  = C_5_R*c12dx1
  c2dx1   = C_2_R   *vmult_imin

  c12dxmx = C_05_R  *vmult_imax
  c52dxmx = C_5_R*c12dxmx
  c2dxmx  = C_2_R   *vmult_imax

!  c0304 = C_075_R
  c0118 = C_1_R/C_18_R
  c0102 = C_1_R/C_2_R
  c0509 = C_5_R/C_9_R

  c1418 = C_14_R/C_18_R
  c0136 = C_1_R /C_36_R

  DO jk =srt,end
! third-order
     d(jk,1)    =-c52dx1 *u(jk,1)    + c2dx1 *u(jk,2)      + c12dx1 *u(jk,3)
     d(jk,imax) = c52dxmx*u(jk,imax) - c2dxmx*u(jk,imax-1) - c12dxmx*u(jk,imax-2)
! fourth-order centered
!     d(jk,2)      = c0304 *(u(jk,3)    - u(jk,1)     )
!     d(jk,imax-1) = c0304 *(u(jk,imax) - u(jk,imax-2))
! fifth-order biased
     d(jk,2)      = c0118*u(jk,4)      +u(jk,3)      -c0102*u(jk,2)      -c0509*u(jk,1)
     d(jk,imax-1) =-c0118*u(jk,imax-3) -u(jk,imax-2) +c0102*u(jk,imax-1) +c0509*u(jk,imax)
  ENDDO

  DO i = 3,imax-2
     DO jk = srt,end
        d(jk,i) = c1418*(u(jk,i+1) - u(jk,i-1)) + c0136*(u(jk,i+2) - u(jk,i-2))
     ENDDO
  ENDDO
!$omp end parallel 
  RETURN
END SUBROUTINE FDM_C1N6_RHS

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2010/11/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Same as before, but imposing a given BC
!#
!########################################################################
!# ARGUMENTS 
!#
!# jkmax In    number of systems to solve
!# ibc   In    BCs: 1 u'_1 given
!#                  2 u'_N given
!#
!# d     Out   right-hand side vector of the linear system
!#
!########################################################################

!########################################################################
!Left-hand side; tridiagonal matrix of the linear system
!########################################################################
SUBROUTINE FDM_C1N6_BCS_LHS(imax, ibc, dx, a,b,c)

  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax, ibc
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i
  TREAL c0103

! #######################################################################
  c0103 = C_1_R/C_3_R

! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
  a(1     ) = C_0_R ! padding
  b(1     ) = C_1_R *dx(1)
  c(1     ) = C_2_R *dx(2)
! fourth-order centered
!   a(2     ) = C_025_R *dx(1)
!   b(2     ) = C_1_R   *dx(2)
!   c(2     ) = C_025_R *dx(3)
! fifth-order biased
  a(2) = C_1_R/C_6_R *dx(1)
  b(2) = C_1_R       *dx(2)
  c(2) = C_1_R/C_2_R *dx(3)
  DO i = 3,imax-2
  a(i     ) = c0103 *dx(i-1)
  b(i     ) = C_1_R *dx(i  )
  c(i     ) = c0103 *dx(i+1)
  ENDDO
! fourth-order centered
!   a(imax-1) = C_025_R *dx(imax-2)
!   b(imax-1) = C_1_R   *dx(imax-1)
!   c(imax-1) = C_025_R *dx(imax  )
!fifth-order biased
  a(imax-1) = C_1_R/C_2_R *dx(imax-2)
  b(imax-1) = C_1_R       *dx(imax-1)
  c(imax-1) = C_1_R/C_6_R *dx(imax  )
! third-order biased
  a(imax  ) = C_2_R *dx(imax-1)
  b(imax  ) = C_1_R *dx(imax  )
  c(imax  ) = C_0_R ! padding

! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
  IF ( ibc .EQ. 1 .OR. ibc .EQ. 3 ) THEN
     b(1) =-b(1)*C_2_R/C_5_R            ! first row
     c(1) =-c(1)*C_2_R/C_5_R
!      a(1) = C_1_R/(C_10_R*C_2_R) *dx(1) !-a21R fourth-order centered
!      a(2) = C_1_R/(C_10_R*C_9_R) *dx(1)
     a(1) = C_1_R/(C_9_R *C_2_R) *dx(1) !-a21R fifth-order biased
     a(2) = C_1_R/(C_10_R*C_9_R) *dx(1)
!      b(2) = C_2_R / C_5_R        *dx(2) ! A22R fourth-order centered
!      a(3) = C_14_R/(C_5_R*C_9_R) *dx(2)
     b(2) = C_5_R / C_9_R        *dx(2) ! A22R fifth-order biased
     a(3) = C_14_R/(C_5_R*C_9_R) *dx(2)

  ENDIF
  IF ( ibc .EQ. 2 .OR. ibc .EQ. 3 ) THEN
!      c(imax-2) = C_14_R/(C_5_R*C_9_R) *dx(imax-1) ! A11R fourth-order centered
!      b(imax-1) = C_2_R / C_5_R        *dx(imax-1)
     c(imax-2) = C_14_R/(C_5_R*C_9_R) *dx(imax-1) ! A11R fifth-order biased
     b(imax-1) = C_5_R / C_9_R        *dx(imax-1)
!      c(imax  ) = C_1_R/(C_10_R*C_2_R) *dx(imax)   !-a1nR fourth-order centered
!      c(imax-1) = C_1_R/(C_10_R*C_9_R) *dx(imax)
     c(imax  ) = C_1_R/(C_9_R *C_2_R) *dx(imax)   !-a1nR fifth-order biased
     c(imax-1) = C_1_R/(C_10_R*C_9_R) *dx(imax)
     b(imax  ) = b(imax)*C_2_R/C_5_R              ! last row
     a(imax  ) = a(imax)*C_2_R/C_5_R

  ENDIF

  RETURN
END SUBROUTINE FDM_C1N6_BCS_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N6_BCS_RHS(imax,jkmax, ibc, u,d)

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax,jkmax, ibc
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, imin_loc, imax_loc
  TREAL c0105, c0102, c0305, c0304, c0405, c0502, c0180, c0140, c0144
  TREAL c0709, c0136
  TREAL c0509, c0809, c0118, c1718

! #######################################################################
  c0105 = C_1_R/C_5_R
  c0102 = C_1_R/C_2_R
  c0305 = C_3_R/C_5_R
  c0304 = C_3_R/C_4_R
  c0405 = C_4_R/C_5_R
  c0502 = C_5_R/C_2_R
  c0180 = C_1_R/(C_10_R*C_18_R)
  c0140 = C_14_R*C_10_R
  c0144 = C_12_R*C_12_R

  c0709 = C_7_R/C_9_R
  c0136 = C_1_R/C_36_R
  c0509 = C_5_R /C_9_R
  c0809 = C_8_R /C_9_R
  c0118 = C_1_R /C_18_R
  c1718 = C_17_R/C_18_R

  imin_loc = 3; imax_loc = imax-2

! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
  IF      ( ibc .EQ. 1 ) THEN
! contribution to u_1
     d(:,1     ) = c0405*u(:,2) +c0105*u(:,3)

! B22R
!      d(:,2     ) = c0305*(u(:,3)-u(:,2)) !4th
!      d(:,3     ) = c0180*(C_5_R*u(:,5)+c0140*u(:,4)-u(:,3)-c0144*u(:,2))
     d(:,2     ) = c0118*u(:,4) +c0809*u(:,3) -c1718*u(:,2) !5th
     d(:,3     ) = c0180*(C_5_R*u(:,5)+c0140*u(:,4)-u(:,3)-c0144*u(:,2))
     imin_loc = imin_loc + 1

     d(:,imax-1) =-c0118*u(:,imax-3) -u(:,imax-2) +c0102*u(:,imax-1) +c0509*u(:,imax) !5th
!     d(:,imax-1) = c0304*(u(:,imax)-u(:,imax-2))                         !4th
     d(:,imax  ) = c0502*u(:,imax) -C_2_R*u(:,imax-1) -c0102*u(:,imax-2) !3rd

! -------------------------------------------------------------------
  ELSE IF ( ibc .EQ. 2 ) THEN
     d(:,1     ) =-c0502*u(:,1) +C_2_R*u(:,2) +c0102*u(:,3) !3rd
!     d(:,2     ) = c0304*(u(:,3)-u(:,1))                    !4th
     d(:,2     ) = c0118*u(:,4) +u(:,3) -c0102*u(:,2) -c0509*u(:,1) !5th

! B11R
     imax_loc = imax_loc - 1
     d(:,imax-2) = c0180*(c0144*u(:,imax-1)+u(:,imax-2)-c0140*u(:,imax-3)-C_5_R*u(:,imax-4))!5th
     d(:,imax-1) = c1718*u(:,imax-1) -c0809*u(:,imax-2) -c0118*u(:,imax-3)
!     d(:,imax-2) = c0180*(c0144*u(:,imax-1)+u(:,imax-2)-c0140*u(:,imax-3)-C_5_R*u(:,imax-4)) !4th
!     d(:,imax-1) = c0305*(u(:,imax-1)-u(:,imax-2))

! contribution to u_n
     d(:,imax  ) = c0405*u(:,imax-1) +c0105*u(:,imax-2)

! -------------------------------------------------------------------
  ELSE IF ( ibc .EQ. 3 ) THEN
! contribution to u_1
     d(:,1     ) = c0405*u(:,2) +c0105*u(:,3)

! B22R
!      d(:,2     ) = c0305*(u(:,3)-u(:,2)) !4th
!      d(:,3     ) = c0180*(C_5_R*u(:,5)+c0140*u(:,4)-u(:,3)-c0144*u(:,2))
     d(:,2     ) = c0118*u(:,4) +c0809*u(:,3) -c1718*u(:,2) !5th
     d(:,3     ) = c0180*(C_5_R*u(:,5)+c0140*u(:,4)-u(:,3)-c0144*u(:,2))
     imin_loc = imin_loc + 1

! B11R
     imax_loc = imax_loc - 1
     d(:,imax-2) = c0180*(c0144*u(:,imax-1)+u(:,imax-2)-c0140*u(:,imax-3)-C_5_R*u(:,imax-4)) !5th
     d(:,imax-1) = c1718*u(:,imax-1) -c0809*u(:,imax-2) -c0118*u(:,imax-3)
!      d(:,imax-2) = c0180*(c0144*u(:,imax-1)+u(:,imax-2)-c0140*u(:,imax-3)-C_5_R*u(:,imax-4)) !4th
!      d(:,imax-1) = c0305*(u(:,imax-1)-u(:,imax-2))

! contribution to u_n
     d(:,imax  ) = c0405*u(:,imax-1) +c0105*u(:,imax-2)

! -------------------------------------------------------------------
  ELSE
     d(:,1     ) =-c0502*u(:,1) +C_2_R*u(:,2) +c0102*u(:,3) !3rd
!     d(:,2     ) = c0304*(u(:,3)-u(:,1))                    !4th
     d(:,2     ) = c0118*u(:,4) +u(:,3) -c0102*u(:,2) -c0509*u(:,1) !5th

     d(:,imax-1) =-c0118*u(:,imax-3) -u(:,imax-2) +c0102*u(:,imax-1) +c0509*u(:,imax) !5th
!     d(:,imax-1) = c0304*(u(:,imax)-u(:,imax-2))                         !4th
     d(:,imax  ) = c0502*u(:,imax) -C_2_R*u(:,imax-1) -c0102*u(:,imax-2) !3rd

  ENDIF

! -------------------------------------------------------------------
! Interior points
! -------------------------------------------------------------------
  DO i = imin_loc,imax_loc
     d(:,i) = c0709*(u(:,i+1) - u(:,i-1)) + c0136*(u(:,i+2) - u(:,i-2))
  ENDDO
 
  RETURN
END SUBROUTINE FDM_C1N6_BCS_RHS

