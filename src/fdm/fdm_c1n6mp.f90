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
!# 6th order pentadiagonal compact scheme by JCP Lele 1992, periodic.
!# Interior points according to Eq. 2.1.10. Similar truncation error like 
!# Eq. 2.1.7 with (\alpha=1/3). Here the alpha value is chosen such, 
!# that no inflection point in w'(w) appears.
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
#define C1N6M_ALPHA 0.604730585697398d+0
#define C1N6M_BETA  0.108558900945626d+0
! coefficients RHS
#define C1N6M_AD2   0.619462713898740d+0
#define C1N6M_BD4   0.284700510015759d+0
#define C1N6M_CD6   0.814191757092195d-2

! #######################################################################
! Left-hand side; pentadiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C1N6MP_LHS(imax, dx, a,b,c,d,e)

  ! USE TLAB_VARS, ONLY : C1N6M_ALPHA, C1N6M_BETA
  
  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c,d,e

! -------------------------------------------------------------------
  TINTEGER                            :: i

! ###################################################################
  DO i = 1,imax
    a(i) = C1N6M_BETA 
    b(i) = C1N6M_ALPHA
    c(i) = C_1_R
    d(i) = C1N6M_ALPHA
    e(i) = C1N6M_BETA 
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
END SUBROUTINE FDM_C1N6MP_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N6MP_RHS(imax,jkmax, u,d)

  ! USE TLAB_VARS, ONLY : C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
  
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                  :: i, jk
  TINTEGER                                  :: im3, im2, im1, ip1, ip2, ip3

! #######################################################################
  DO i = 1, imax
    im3 = MOD(i+imax-4, imax) + 1 ! n-3 
    im2 = MOD(i+imax-3, imax) + 1 ! n-2 
    im1 = MOD(i+imax-2, imax) + 1 ! n-1 
    ip1 = MOD(i,        imax) + 1 ! n+1 
    ip2 = MOD(i+1,      imax) + 1 ! n+2 
    ip3 = MOD(i+2,      imax) + 1 ! n+3 
    DO jk = 1,jkmax
      d(jk,i) = C1N6M_AD2 * (u(jk,ip1) - u(jk,im1)) + &
                C1N6M_BD4 * (u(jk,ip2) - u(jk,im2)) + &
                C1N6M_CD6 * (u(jk,ip3) - u(jk,im3))  
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE FDM_C1N6MP_RHS