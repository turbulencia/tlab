#include "types.h"

#define C_12D11_L .1090909090909090d+1
#define C_51D22_L .2318181818181818d+1
#define C_03D44_L .6818181818181818d-1
#define C_02D11_L .1818181818181818d+0

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2012/12/20 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic,
!# in order to solve the BVP
!#
!#     u''_i - \lamba^2 u_i = f_i  N-2 eqns
!#     u_1 and u_N given           2   eqns
!#     Au'' = Bu                   N   eqns
!#
!# The system of N-2 eqns:
!# 
!#                    (B - \lambda^2 A)u = Af = g
!#
!# is established in this routine, giving diagonals a-e and g (see notes).
!# Interior points 6th-order according to Eq. 2.2.7.
!# The second point from Eq. 2.2.6 forth-order (b=0).
!# The first point from third-order biased Eq. 4.3.1.
!#
!# Solution array does not appear in this routine.
!#
!########################################################################

!########################################################################
!Left-hand side; pentadiagonal matrix of the linear system and f1 and f2
!########################################################################
SUBROUTINE INT_C2N6_LHS_E(imax, dx, lambda2, a,b,c,d,e, f1,f2)

  IMPLICIT NONE

  TREAL lambda2
  TINTEGER,                 INTENT(IN)  :: imax       ! original size; here using only 2:imax-1
  TREAL, DIMENSION(imax,2), INTENT(IN)  :: dx
  TREAL, DIMENSION(imax),   INTENT(OUT) :: a,b,c,d,e  ! diagonals
  TREAL, DIMENSION(imax),   INTENT(OUT) :: f1,f2      ! forcing term for the hyperbolic sine

! -------------------------------------------------------------------
  TINTEGER i

! ###################################################################
! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
  a(2     ) = C_0_R ! padding
  b(2     ) = C_0_R ! padding
  c(2     ) =-C_3_R - lambda2*dx(2,1)*dx(2,1)
  d(2     ) = C_3_R + lambda2*dx(3,1)*dx(3,1)
  e(2     ) =-C_1_R      
! fourth-order centered
  a(3     ) = C_0_R ! padding
  b(3     ) = C_12D11_L - C_02D11_L*lambda2*dx(2,1)*dx(2,1)
  c(3     ) =-C_51D22_L -           lambda2*dx(3,1)*dx(3,1)
  d(3     ) = C_12D11_L - C_02D11_L*lambda2*dx(4,1)*dx(4,1)   
  e(3     ) = C_03D44_L
! sixth-order centered
  DO i = 4,imax-3
  a(i     ) = C_03D44_L
  b(i     ) = C_12D11_L - C_02D11_L*lambda2*dx(i-1,1)*dx(i-1,1)
  c(i     ) =-C_51D22_L -           lambda2*dx(i  ,1)*dx(i  ,1)
  d(i     ) = C_12D11_L - C_02D11_L*lambda2*dx(i+1,1)*dx(i+1,1)
  e(i     ) = C_03D44_L
  ENDDO
! fourth-order centered
  a(imax-2) = C_03D44_L
  b(imax-2) = C_12D11_L - C_02D11_L*lambda2*dx(imax-3,1)*dx(imax-3,1)
  c(imax-2) =-C_51D22_L -           lambda2*dx(imax-2,1)*dx(imax-2,1)
  d(imax-2) = C_12D11_L - C_02D11_L*lambda2*dx(imax-1,1)*dx(imax-1,1)
  e(imax-2) = C_0_R ! padding
! third-order biased
  a(imax-1) =-C_1_R
  b(imax-1) = C_3_R + lambda2*dx(imax-2,1)*dx(imax-2,1)
  c(imax-1) =-C_3_R - lambda2*dx(imax-1,1)*dx(imax-1,1)
  d(imax-1) = C_0_R ! padding
  e(imax-1) = C_0_R ! padding
    
! -------------------------------------------------------------------
! Setting the RHS for the hyperbolic sine 
! The minus sign in included here to save ops
! -------------------------------------------------------------------
  f1 = C_0_R
  f1(1     ) = C_1_R ! This element is simply the solution at imin of s(-)
  f1(2     ) =-C_1_R
  f1(3     ) =-C_03D44_L

  f2 = C_0_R
  f2(imax-2) =-C_03D44_L
  f2(imax-1) =-C_1_R
  f2(imax  ) = C_1_R ! This element is simply the solution at imax of s(+)

  RETURN
END SUBROUTINE INT_C2N6_LHS_E

! #######################################################################
! Right-hand side; jkmax forcing terms at the same time
! #######################################################################
SUBROUTINE INT_C2N6_RHS(imax,jkmax, dx, f,g)

  IMPLICIT NONE

  TINTEGER,                    INTENT(IN) :: imax, jkmax
  TREAL, DIMENSION(imax,2),    INTENT(IN) :: dx
  TREAL, DIMENSION(jkmax,imax),INTENT(IN) :: f
  TREAL, DIMENSION(jkmax,imax),INTENT(OUT):: g

! -------------------------------------------------------------------
  TINTEGER i

! ###################################################################
! Boundary conditions
  g(:,1     ) = C_0_R ! This element is simply the solution at imin of p(0)
  g(:,2     ) = f(:,2     )*dx(2,1)     *dx(2,1)      - f(:,3)     *dx(3,1)     *dx(3,1)
  g(:,imax-1) = f(:,imax-1)*dx(imax-1,1)*dx(imax-1,1) - f(:,imax-2)*dx(imax-2,1)*dx(imax-2,1)
  g(:,imax  ) = C_0_R ! This element is simply the solution at imax of p(0)

! Interior points
  DO i = 3,imax-2
  g(:,i     ) = f(:,i)*dx(i,1)*dx(i,1) + &
       C_02D11_L*( f(:,i-1)*dx(i-1,1)*dx(i-1,1) + f(:,i+1)*dx(i+1,1)*dx(i+1,1) )
  ENDDO

  RETURN
END SUBROUTINE INT_C2N6_RHS
