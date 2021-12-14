!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/13 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!# 
!# Implementation of the first derivative finite difference with
!# 6th order pentadiagonal compact scheme by JCP Lele 1992, non-periodic.
!# Interior points according to Eq. 2.1.10. Similar truncation error like 
!# Eq. 2.1.7 with (\alpha=1/3). Here alpha value is chosen such, 
!# that no inflection point in w'(w) appears.
!# Carpenter, Gottlieb and Aberbanel, JCP, 1993, study the effect of 
!# boundary points on stability. Scheme 3-5-c6--penta c6--c6-5-3 
!# is implmented.
!#
!########################################################################
!# ARGUMENTS 
!#
!# u    In    function to be diferentiated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################
#include "types.h"

! coefficients LHS
#define C_C1N6MP_ALPHA_L 0.604730585697398d+0
#define C_C1N6MP_BETA_L  0.108558900945626d+0
! coefficients RHS
#define C_C1N6MP_AD2_L   0.619462713898740d+0
#define C_C1N6MP_BD4_L   0.284700510015759d+0
#define C_C1N6MP_CD6_L   0.814191757092195d-2

!########################################################################
! Left-hand side; pentadiagonal matrix of the linear system
!########################################################################
SUBROUTINE FDM_C1N6M_LHS(imax, imin_set_zero,imax_set_zero, dx,a,b,c,d,e)

  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax, imin_set_zero,imax_set_zero
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c,d,e

! -------------------------------------------------------------------
  TINTEGER                            :: i
  TREAL                               :: vmult_imin, vmult_imax, c0103

! #######################################################################
  c0103 = C_1_R/C_3_R

  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! #######################################################################

! third-order biased
  a(1)    = C_0_R
  b(1)    = C_0_R
  c(1)    = C_1_R
  d(1)    = C_2_R * vmult_imin
  e(1)    = C_0_R
  !
  a(imax) = C_0_R
  b(imax) = C_2_R * vmult_imax
  c(imax) = C_1_R
  d(imax) = C_0_R
  e(imax) = C_0_R

! fifth-order biased
  a(2)      = C_0_R
  b(2)      = C_1_R/C_6_R
  c(2)      = C_1_R
  d(2)      = C_1_R/C_2_R
  e(2)      = C_0_R
  !
  a(imax-1) = C_0_R
  b(imax-1) = C_1_R/C_2_R
  c(imax-1) = C_1_R
  d(imax-1) = C_1_R/C_6_R
  e(imax-1) = C_0_R
  
! sixth-order centered (alpha=1/3)
  a(3)      = C_0_R
  b(3)      = c0103
  c(3)      = C_1_R
  d(3)      = c0103
  e(3)      = C_0_R
  !
  a(imax-2) = C_0_R
  b(imax-2) = c0103
  c(imax-2) = C_1_R
  d(imax-2) = c0103
  e(imax-2) = C_0_R
  
! sixth-order modified centered
  DO i = 4,imax-3
    a(i) = C_C1N6MP_BETA_L 
    b(i) = C_C1N6MP_ALPHA_L
    c(i) = C_1_R
    d(i) = C_C1N6MP_ALPHA_L
    e(i) = C_C1N6MP_BETA_L 
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
  e(imax-1) = e(imax-1) * dx(1)
  d(imax  ) = d(imax  ) * dx(1)
  c(1     ) = c(1     ) * dx(1)
  b(2     ) = b(2     ) * dx(1)
  a(3     ) = a(3     ) * dx(1)
  !
  e(imax  ) = e(imax  ) * dx(2)
  d(1     ) = d(1     ) * dx(2)
  c(2     ) = c(2     ) * dx(2)
  b(3     ) = b(3     ) * dx(2)
  a(4     ) = a(4     ) * dx(2)
  !
  Do i = 3, imax-2
    e(i-2) = e(i-2) * dx(i)
    d(i-1) = d(i-1) * dx(i)
    c(i  ) = c(i  ) * dx(i)
    b(i+1) = b(i+1) * dx(i)
    a(i+2) = a(i+2) * dx(i)
  ENDDO
  !
  e(imax-3) = e(imax-3) * dx(imax-1)
  d(imax-2) = d(imax-2) * dx(imax-1)
  c(imax-1) = c(imax-1) * dx(imax-1)
  b(imax  ) = b(imax  ) * dx(imax-1)
  a(1     ) = a(1     ) * dx(imax-1)
  !
  e(imax-2) = e(imax-2) * dx(imax)
  d(imax-1) = d(imax-1) * dx(imax)
  c(imax  ) = c(imax  ) * dx(imax)
  b(1     ) = b(1     ) * dx(imax)
  a(2     ) = a(2     ) * dx(imax)

  RETURN
END SUBROUTINE FDM_C1N6M_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N6M_RHS(imax,jkmax, imin_set_zero,imax_set_zero, u,d)

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax,jkmax, imin_set_zero,imax_set_zero
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                  :: i, jk
  TREAL                                     :: vmult_imin, vmult_imax
  TREAL                                     :: c1418, c0136, c0118, c0102, c0509
  TREAL                                     :: c52dx1, c2dx1, c12dx1, c52dxmx, c2dxmx, c12dxmx

! #######################################################################
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_set_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! #######################################################################

  c12dx1  = C_05_R  *vmult_imin
  c52dx1  = C_5_R   *c12dx1
  c2dx1   = C_2_R   *vmult_imin

  c12dxmx = C_05_R  *vmult_imax
  c52dxmx = C_5_R   *c12dxmx
  c2dxmx  = C_2_R   *vmult_imax

  c0118   = C_1_R/C_18_R
  c0102   = C_1_R/C_2_R
  c0509   = C_5_R/C_9_R

  c1418   = C_14_R/C_18_R
  c0136   = C_1_R /C_36_R


  DO jk = 1,jkmax
    ! third-order
    d(jk,1)      =-c52dx1 *u(jk,1)    + c2dx1 *u(jk,2)      + c12dx1 *u(jk,3)
    d(jk,imax)   = c52dxmx*u(jk,imax) - c2dxmx*u(jk,imax-1) - c12dxmx*u(jk,imax-2)
    ! fifth-order biased
    d(jk,2)      = c0118*u(jk,4)      +u(jk,3)      -c0102*u(jk,2)      -c0509*u(jk,1)
    d(jk,imax-1) =-c0118*u(jk,imax-3) -u(jk,imax-2) +c0102*u(jk,imax-1) +c0509*u(jk,imax)
    ! 6th-order centered with alpha=(1/3)
    d(jk,3)      = c1418*(u(jk,     4) - u(jk,     2)) + c0136*(u(jk,    5) - u(jk,     1))
    d(jk,imax-2) = c1418*(u(jk,imax-1) - u(jk,imax-3)) + c0136*(u(jk, imax) - u(jk,imax-4))
  ENDDO
    
  DO i = 4,imax-3
    DO jk = 1,jkmax
      d(jk,i) = C_C1N6MP_AD2_L * (u(jk,i+1) - u(jk,i-1)) + &
                C_C1N6MP_BD4_L * (u(jk,i+2) - u(jk,i-2)) + &
                C_C1N6MP_CD6_L * (u(jk,i+3) - u(jk,i-3))  
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE FDM_C1N6M_RHS

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/14 - J. Kostelecky
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
SUBROUTINE FDM_C1N6M_BCS_LHS(imax, ibc, dx, a,b,c,d,e)

  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax, ibc
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c,d,e

! -------------------------------------------------------------------
  TINTEGER                            :: i
  TREAL                               :: c0103

! #######################################################################
  c0103 = C_1_R/C_3_R

  ! third-order biased
  a(1)    = C_0_R
  b(1)    = C_0_R
  c(1)    = C_1_R
  d(1)    = C_2_R
  e(1)    = C_0_R
  !
  a(imax) = C_0_R
  b(imax) = C_2_R
  c(imax) = C_1_R
  d(imax) = C_0_R
  e(imax) = C_0_R

! fifth-order biased
  a(2)      = C_0_R
  b(2)      = C_1_R/C_6_R
  c(2)      = C_1_R
  d(2)      = C_1_R/C_2_R
  e(2)      = C_0_R
  !
  a(imax-1) = C_0_R
  b(imax-1) = C_1_R/C_2_R
  c(imax-1) = C_1_R
  d(imax-1) = C_1_R/C_6_R
  e(imax-1) = C_0_R
  
! sixth-order centered (alpha=1/3)
  a(3)      = C_0_R
  b(3)      = c0103
  c(3)      = C_1_R
  d(3)      = c0103
  e(3)      = C_0_R
  !
  a(imax-2) = C_0_R
  b(imax-2) = c0103
  c(imax-2) = C_1_R
  d(imax-2) = c0103
  e(imax-2) = C_0_R
  
! sixth-order modified centered
  DO i = 4,imax-3
    a(i) = C_C1N6MP_BETA_L 
    b(i) = C_C1N6MP_ALPHA_L
    c(i) = C_1_R
    d(i) = C_C1N6MP_ALPHA_L
    e(i) = C_C1N6MP_BETA_L 
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
  e(imax-1) = e(imax-1) * dx(1)
  d(imax  ) = d(imax  ) * dx(1)
  c(1     ) = c(1     ) * dx(1)
  b(2     ) = b(2     ) * dx(1)
  a(3     ) = a(3     ) * dx(1)
  !
  e(imax  ) = e(imax  ) * dx(2)
  d(1     ) = d(1     ) * dx(2)
  c(2     ) = c(2     ) * dx(2)
  b(3     ) = b(3     ) * dx(2)
  a(4     ) = a(4     ) * dx(2)
  !
  Do i = 3, imax-2
    e(i-2) = e(i-2) * dx(i)
    d(i-1) = d(i-1) * dx(i)
    c(i  ) = c(i  ) * dx(i)
    b(i+1) = b(i+1) * dx(i)
    a(i+2) = a(i+2) * dx(i)
  ENDDO
  !
  e(imax-3) = e(imax-3) * dx(imax-1)
  d(imax-2) = d(imax-2) * dx(imax-1)
  c(imax-1) = c(imax-1) * dx(imax-1)
  b(imax  ) = b(imax  ) * dx(imax-1)
  a(1     ) = a(1     ) * dx(imax-1)
  !
  e(imax-2) = e(imax-2) * dx(imax)
  d(imax-1) = d(imax-1) * dx(imax)
  c(imax  ) = c(imax  ) * dx(imax)
  b(1     ) = b(1     ) * dx(imax)
  a(2     ) = a(2     ) * dx(imax)

! -------------------------------------------------------------------
! Boundary conditions
! -------------------------------------------------------------------
  IF ( ibc .EQ. 1 .OR. ibc .EQ. 3 ) THEN
    ! first row
    c(1) = - c(1) * C_2_R/C_5_R            
    d(1) = - d(1) * C_2_R/C_5_R
    ! A21R fifth-order biased
    b(1) = C_1_R  / (C_9_R *C_2_R) * dx(1)
    b(2) = C_1_R  / (C_10_R*C_9_R) * dx(1)
    ! A22R fifth-order biased
    c(2) = C_5_R  /  C_9_R         * dx(2) 
    b(3) = C_14_R / (C_5_R*C_9_R)  * dx(2)
  ENDIF

  IF ( ibc .EQ. 2 .OR. ibc .EQ. 3 ) THEN
    ! A11R fifth-order biased
    d(imax-2) = C_14_R / (C_5_R*C_9_R)  * dx(imax-1) 
    c(imax-1) = C_5_R  /  C_9_R         * dx(imax-1)
    ! A1nR fifth-order biased 
    d(imax  ) = C_1_R  / (C_9_R *C_2_R) * dx(imax)   
    d(imax-1) = C_1_R  / (C_10_R*C_9_R) * dx(imax)
    ! last row
    c(imax  ) = c(imax) * C_2_R/C_5_R              
    b(imax  ) = b(imax) * C_2_R/C_5_R
  ENDIF

  RETURN
END SUBROUTINE FDM_C1N6M_BCS_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N6M_BCS_RHS(imax,jkmax, ibc, u,d)

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax,jkmax, ibc
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER :: i, imin_loc, imax_loc
  TREAL    :: c0105, c0102, c0305, c0304, c0405, c0502, c0180, c0140, c0144
  TREAL    :: c0709, c0136
  TREAL    :: c0509, c0809, c0118, c1718

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
    d(:,2     ) = c0118*u(:,4) +c0809*u(:,3) -c1718*u(:,2) !5th
    d(:,3     ) = c0180*(C_5_R*u(:,5)+c0140*u(:,4)-u(:,3)-c0144*u(:,2))

    imin_loc = imin_loc + 1

    d(:,imax-1) =-c0118*u(:,imax-3) -u(:,imax-2) +c0102*u(:,imax-1) +c0509*u(:,imax) !5th
    d(:,imax  ) = c0502*u(:,imax) -C_2_R*u(:,imax-1) -c0102*u(:,imax-2) !3rd
! -------------------------------------------------------------------
  ELSE IF ( ibc .EQ. 2 ) THEN
    d(:,1     ) =-c0502*u(:,1) +C_2_R*u(:,2) +c0102*u(:,3) !3rd
    d(:,2     ) = c0118*u(:,4) +u(:,3) -c0102*u(:,2) -c0509*u(:,1) !5th

    ! B11R
    imax_loc = imax_loc - 1
    d(:,imax-2) = c0180*(c0144*u(:,imax-1)+u(:,imax-2)-c0140*u(:,imax-3)-C_5_R*u(:,imax-4))!5th
    d(:,imax-1) = c1718*u(:,imax-1) -c0809*u(:,imax-2) -c0118*u(:,imax-3)

    ! contribution to u_n
    d(:,imax  ) = c0405*u(:,imax-1) +c0105*u(:,imax-2)
! -------------------------------------------------------------------
  ELSE IF ( ibc .EQ. 3 ) THEN
    ! contribution to u_1
    d(:,1     ) = c0405*u(:,2) +c0105*u(:,3)

    ! B22R
    d(:,2     ) = c0118*u(:,4) +c0809*u(:,3) -c1718*u(:,2) !5th
    d(:,3     ) = c0180*(C_5_R*u(:,5)+c0140*u(:,4)-u(:,3)-c0144*u(:,2))
    imin_loc = imin_loc + 1

    ! B11R
    imax_loc = imax_loc - 1
    d(:,imax-2) = c0180*(c0144*u(:,imax-1)+u(:,imax-2)-c0140*u(:,imax-3)-C_5_R*u(:,imax-4)) !5th
    d(:,imax-1) = c1718*u(:,imax-1) -c0809*u(:,imax-2) -c0118*u(:,imax-3)

    ! contribution to u_n
    d(:,imax  ) = c0405*u(:,imax-1) +c0105*u(:,imax-2)
! -------------------------------------------------------------------
  ELSE
    d(:,1     ) =-c0502*u(:,1) +C_2_R*u(:,2) +c0102*u(:,3) !3rd
    d(:,2     ) = c0118*u(:,4) +u(:,3) -c0102*u(:,2) -c0509*u(:,1) !5th

    d(:,imax-1) =-c0118*u(:,imax-3) -u(:,imax-2) +c0102*u(:,imax-1) +c0509*u(:,imax) !5th
    d(:,imax  ) = c0502*u(:,imax) -C_2_R*u(:,imax-1) -c0102*u(:,imax-2) !3rd
 
  ENDIF

! -------------------------------------------------------------------
! Interior points
! -------------------------------------------------------------------
  ! 6th-order centered with alpha=(1/3)
  i = imin_loc
  d(:,i) = c0709*(u(:,i+1) - u(:,i-1)) + c0136*(u(:,i+2) - u(:,i-2))
  i = imax_loc
  d(:,i) = c0709*(u(:,i+1) - u(:,i-1)) + c0136*(u(:,i+2) - u(:,i-2))
  ! 6th-order centered - pentadiagonal
  DO i = imin_loc+1, imax_loc-1
    d(:,i) = C_C1N6MP_AD2_L * (u(:,i+1) - u(:,i-1)) + &
             C_C1N6MP_BD4_L * (u(:,i+2) - u(:,i-2)) + &
             C_C1N6MP_CD6_L * (u(:,i+3) - u(:,i-3))  
  ENDDO

  RETURN
END SUBROUTINE FDM_C1N6M_BCS_RHS