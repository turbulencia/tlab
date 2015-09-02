#include "types.h"

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2013/01/25 - J.P. Mellado
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
SUBROUTINE INT_C2N6N_LHS_E(imax, lhs,rhs, lambda2, a,b,c,d,e, f1,f2)

  IMPLICIT NONE

  TREAL lambda2
  TINTEGER,                 INTENT(IN)  :: imax       ! original size; here using only 2:imax-1
  TREAL, DIMENSION(imax,3), INTENT(IN)  :: lhs
  TREAL, DIMENSION(imax,4), INTENT(IN)  :: rhs
  TREAL, DIMENSION(imax),   INTENT(OUT) :: a,b,c,d,e  ! diagonals
  TREAL, DIMENSION(imax),   INTENT(OUT) :: f1,f2      ! forcing term for the hyperbolic sine

! -------------------------------------------------------------------
  TINTEGER i
  TREAL dummy1, dummy2

! ###################################################################
! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system (array C22R)
! -------------------------------------------------------------------
  i = 2; dummy1 = lhs(2,1)/lhs(1,2)
  a(i) = C_0_R ! padding
  b(i) = C_0_R ! padding
  c(i) = C_1_R    - lambda2 *lhs(i,2) - dummy1 *( rhs(1,3) - lambda2 *lhs(1,3) )
  d(i) = rhs(i,3) - lambda2 *lhs(i,3) - dummy1 *  rhs(1,4) 
  e(i) = rhs(i,4)                     - dummy1 *  rhs(1,1)

  i = 3
  a(i) = C_0_R ! padding
  b(i) = rhs(i,2) - lambda2 *lhs(i,1)
  c(i) = C_1_R    - lambda2 *lhs(i,2)
  d(i) = rhs(i,3) - lambda2 *lhs(i,3)
  e(i) = rhs(i,4)

  DO i = 4,imax-3
  a(i) = rhs(i,1)
  b(i) = rhs(i,2) - lambda2 *lhs(i,1)
  c(i) = C_1_R    - lambda2 *lhs(i,2)
  d(i) = rhs(i,3) - lambda2 *lhs(i,3)
  e(i) = rhs(i,4)
  ENDDO

  i = imax-2
  a(i) = rhs(i,1)
  b(i) = rhs(i,2) - lambda2 *lhs(i,1)
  c(i) = C_1_R    - lambda2 *lhs(i,2)
  d(i) = rhs(i,3) - lambda2 *lhs(i,3)
  e(i) = C_0_R

  i = imax-1; dummy2 = lhs(imax-1,3)/lhs(imax,2)
  a(i) = rhs(i,1)                     - dummy2 *  rhs(imax,4)
  b(i) = rhs(i,2) - lambda2 *lhs(i,1) - dummy2 *  rhs(imax,1)
  c(i) = C_1_R    - lambda2 *lhs(i,2) - dummy2 *( rhs(imax,2) - lambda2 *lhs(imax,1) )
  d(i) = C_0_R ! padding
  e(i) = C_0_R ! padding
    
! -------------------------------------------------------------------
! Setting the RHS for the hyperbolic sine 
! The minus sign in included here to save ops
! -------------------------------------------------------------------
  f1 = C_0_R ! b21R
  f1(1     ) = C_1_R    ! This element is simply the solution at imin of s(-)
  f1(2     ) =-(rhs(2,2)-dummy1)
  f1(3     ) =- rhs(3,1)

  f2 = C_0_R ! b2nR
  f2(imax-2) =- rhs(imax-2,4)
  f2(imax-1) =-(rhs(imax-1,3)-dummy2)
  f2(imax  ) = C_1_R ! This element is simply the solution at imax of s(+)

  RETURN
END SUBROUTINE INT_C2N6N_LHS_E

! #######################################################################
! Right-hand side; mmax forcing terms at the same time
! #######################################################################
SUBROUTINE INT_C2N6N_RHS(imax,mmax, lhs, f,g)

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN) :: imax, mmax
  TREAL, DIMENSION(imax,3),   INTENT(IN) :: lhs
  TREAL, DIMENSION(mmax,imax),INTENT(IN) :: f
  TREAL, DIMENSION(mmax,imax),INTENT(OUT):: g

! -------------------------------------------------------------------
  TINTEGER i
  TREAL dummy

! ###################################################################
  g(:,1) = C_0_R ! This element is simply the solution at imin of p(0)

  i = 2;      dummy = lhs(i,2) - lhs(2,1)/lhs(1,2)*lhs(1,3)
  g(:,i) =                     f(:,i  )*dummy    + f(:,i+1)*lhs(i,3)
  DO i = 3,imax-2 ! Interior points
  g(:,i) = f(:,i-1)*lhs(i,1) + f(:,i  )*lhs(i,2) + f(:,i+1)*lhs(i,3)
  ENDDO
  i = imax-1; dummy = lhs(i,2) - lhs(imax-1,3)/lhs(imax,2)*lhs(imax,1)
  g(:,i) = f(:,i-1)*lhs(i,1) + f(:,i  )*dummy

  g(:,imax  ) = C_0_R ! This element is simply the solution at imax of p(0)

  RETURN
END SUBROUTINE INT_C2N6N_RHS
