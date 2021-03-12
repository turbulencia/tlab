#include "types.h"

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2010/10/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic,
!# in order to solve the IVP
!#
!#     u'_i + \lamba u_i = f_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!# The system of N-1 eqns:
!# 
!#                    (B + \lambda A)u = Af = g
!#
!# is established in this routine, giving diagonals a-e and g (see notes).
!# Interior points 6th-order according to Eq. 2.1.7.
!# The second point from Eq. 2.1.6 forth-order (b=0).
!# The first point from third-order biased Eq. 4.1.3 (d=0).
!#
!# Solution array does not appear in this routine.
!#
!########################################################################
!# ARGUMENTS 
!#
!# jkmax       In    number of systems to solve
!# ibc         In    BCs: 1 u_1 given
!#                        2 u_N given
!#
!# f           Out   forcing term for the exponential
!#
!########################################################################

!########################################################################
!Left-hand side; pentadiagonal matrix of the linear system and f
!########################################################################
SUBROUTINE INT_C1N6_LHS_E(imax, ibc, dx, lambda, a,b,c,d,e, f)

  IMPLICIT NONE

  TREAL lambda
  TINTEGER,               INTENT(IN)  :: imax, ibc
  TREAL, DIMENSION(imax), INTENT(IN)  :: dx
  TREAL, DIMENSION(imax), INTENT(OUT) :: a,b,c,d,e,f

! -------------------------------------------------------------------
  TINTEGER i
  TREAL c0136, c1418, c0103, c0104

! ###################################################################
  c0136 = C_1_R /C_36_R
  c1418 = C_14_R/C_18_R
  c0103 = C_1_R /C_3_R 
  c0104 = C_1_R /C_4_R 

! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
  a(1       ) = C_0_R ! padding
  b(1       ) = C_0_R ! padding
  c(1       ) =-C_5_R/C_2_R + lambda       *dx(1)
  d(1       ) = C_2_R       + lambda*C_2_R *dx(2)
  e(1       ) = C_05_R      
! fourth-order centered
!   a(2       ) = C_0_R ! padding
!   b(2       ) =-C_3_R/C_4_R + lambda*c0104 *dx(1)
!   c(2       ) = C_0_R       + lambda       *dx(2)
!   d(2       ) = C_3_R/C_4_R + lambda*c0104 *dx(3)
!   e(2       ) = C_0_R
! fifth-order biased
  a(2       ) = C_0_R ! padding
  b(2       ) =-C_5_R/C_9_R + lambda*C_1_R/C_6_R *dx(1)
  c(2       ) =-C_1_R/C_2_R + lambda             *dx(2)
  d(2       ) = C_1_R       + lambda*C_1_R/C_2_R *dx(3)
  e(2       ) = C_1_R/C_18_R 
! sixth-order centered
  DO i = 3,imax-2
  a(i       ) =-c0136
  b(i       ) =-c1418       + lambda*c0103  *dx(i-1)
  c(i       ) = C_0_R       + lambda        *dx(i  )
  d(i       ) = c1418       + lambda*c0103  *dx(i+1)
  e(i       ) = c0136
  ENDDO
! fourth-order centered
!   a(imax-1  ) = C_0_R
!   b(imax-1  ) =-C_3_R/C_4_R + lambda*c0104 *dx(imax-2)
!   c(imax-1  ) = C_0_R       + lambda       *dx(imax-1)
!   d(imax-1  ) = C_3_R/C_4_R + lambda*c0104 *dx(imax  )
!   e(imax-1  ) = C_0_R ! padding
! fifth-order biased
  a(  imax-1) =-C_1_R/C_18_R
  b(  imax-1) =-C_1_R       + lambda*C_1_R/C_2_R *dx(imax-2)
  c(  imax-1) = C_1_R/C_2_R + lambda             *dx(imax-1)
  d(  imax-1) = C_5_R/C_9_R + lambda*C_1_R/C_6_R *dx(imax  )
  e(  imax-1) = C_0_R ! padding
! third-order biased
  a(imax    ) =-C_05_R      
  b(imax    ) =-C_2_R       + lambda*C_2_R *dx(imax-1)
  c(imax    ) = C_5_R/C_2_R + lambda       *dx(imax  )
  d(imax    ) = C_0_R ! padding
  e(imax    ) = C_0_R ! padding
    
! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
  IF      ( ibc .EQ. 1 ) THEN
     c(1     ) =-C_5_R/C_2_R                            ! array B22R
!     c(2     ) =-C_1_R/C_2_R + lambda*C_1_R/C_2_R*dx(2) ! array B22R + lambda A22R
!     d(2     ) = C_5_R/C_8_R + lambda*C_1_R/C_4_R*dx(3)
!     e(imax  ) = C_025_R                                 ! element A(imax-1,imax)
     c(2     ) =-C_5_R /C_6_R  + lambda*C_2_R/C_3_R*dx(2) ! B22R + lambda A22R fifth-order biased
     d(2     ) = C_11_R/C_12_R + lambda*C_1_R/C_2_R*dx(3)
     e(imax  ) = C_1_R/C_6_R                              ! element A(imax-1,imax)

  ELSE IF ( ibc .EQ. 2 ) THEN
     a(1)      = C_1_R/C_6_R                                   ! element A(2,1)
     b(imax-1) =-C_11_R/C_12_R + lambda*C_1_R/C_2_R*dx(imax-2) ! B11R + lambda A11R fifth-order biased
     c(imax-1) = C_5_R /C_6_R  + lambda*C_2_R/C_3_R*dx(imax-1)
!     a(1)      = C_025_R                                     ! element A(2,1) 4th
!     b(imax-1) =-C_5_R/C_8_R + lambda*C_1_R/C_4_R*dx(imax-2) ! array B11R + lambda A11R 4th
!     c(imax-1) = C_1_R/C_2_R + lambda*C_1_R/C_2_R*dx(imax-1)
     c(imax  ) = C_5_R/C_2_R                                 ! array B11R

  ENDIF

! -------------------------------------------------------------------
! Setting the RHS for the null space (with minus sign)
! -------------------------------------------------------------------
  f = C_0_R
  IF      ( ibc .EQ. 1 ) THEN
     f(1     ) = C_1_R       ! normalization
!     f(2     ) = C_1_R/C_8_R ! b21R
!     f(3     ) = c0136
     f(2     ) = c0136*C_5_R ! fifth-order biased
     f(3     ) = c0136
  ELSE IF ( ibc .EQ. 2 ) THEN
!     f(imax-2) =-c0136       ! b1NR
!     f(imax-1) =-C_1_R/C_8_R
     f(imax-2) =-c0136       ! fifth-order biased
     f(imax-1) =-c0136*C_5_R
     f(imax  ) = C_1_R       ! normalization
  ENDIF

  RETURN
END SUBROUTINE INT_C1N6_LHS_E

!########################################################################
!Left-hand side case \lambda=0; pentadiagonal matrix of the linear system
!########################################################################
SUBROUTINE INT_C1N6_LHS(imax, ibc, a,b,c,d,e)

  IMPLICIT NONE

  TINTEGER,               INTENT(IN) :: imax, ibc
  TREAL, DIMENSION(imax), INTENT(OUT):: a,b,c,d,e

! -------------------------------------------------------------------
  TREAL c0136, c1418

! ###################################################################
  c0136 = C_1_R /C_36_R
  c1418 = C_14_R/C_18_R

! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
  a(1       ) = C_0_R ! padding
  b(1       ) = C_0_R ! padding
  c(1       ) =-C_5_R/C_2_R
  d(1       ) = C_2_R
  e(1       ) = C_05_R      
! fourth-order centered
!    a(2       ) = C_0_R ! padding
!    b(2       ) =-C_3_R/C_4_R
!    c(2       ) = C_0_R
!    d(2       ) = C_3_R/C_4_R
!    e(2       ) = C_0_R
! fifth-order biased
  a(2       ) = C_0_R ! padding
  b(2       ) =-C_5_R/C_9_R
  c(2       ) =-C_1_R/C_2_R
  d(2       ) = C_1_R
  e(2       ) = C_1_R/C_18_R
! sixth-order centered
  a(3:imax-2) =-c0136
  b(3:imax-2) =-c1418
  c(3:imax-2) = C_0_R
  d(3:imax-2) = c1418
  e(3:imax-2) = c0136
! fourth-order centered
!   a(  imax-1) = C_0_R
!   b(  imax-1) =-C_3_R/C_4_R
!   c(  imax-1) = C_0_R      
!   d(  imax-1) = C_3_R/C_4_R
!   e(  imax-1) = C_0_R ! padding
! fifth-order biased
  a(  imax-1) =-C_1_R/C_18_R
  b(  imax-1) =-C_1_R
  c(  imax-1) = C_1_R/C_2_R      
  d(  imax-1) = C_5_R/C_9_R
  e(  imax-1) = C_0_R ! padding
! third-order biased
  a(  imax  ) =-C_05_R      
  b(  imax  ) =-C_2_R       
  c(  imax  ) = C_5_R/C_2_R 
  d(  imax  ) = C_0_R ! padding
  e(  imax  ) = C_0_R ! padding
     
! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
  IF      ( ibc .EQ. 1 ) THEN ! array B22R
!     c(2     ) =-C_1_R/C_2_R ! fourth-order biased
!     d(2     ) = C_5_R/C_8_R
!     e(imax  ) = C_025_R     ! element A(imax-1,imax)
     c(2     ) =-C_5_R /C_6_R ! fifth-order biased
     d(2     ) = C_11_R/C_12_R
     e(imax  ) = C_1_R/C_6_R
     
  ELSE IF ( ibc .EQ. 2 ) THEN ! array B11R
!     b(imax-1) =-C_5_R/C_8_R ! fourth-order biased
!     c(imax-1) = C_1_R/C_2_R
!     a(1)      = C_025_R     ! element A(2,1)
     b(imax-1) =-C_11_R/C_12_R ! fifth-order biased
     c(imax-1) = C_5_R /C_6_R
     a(1)      = C_1_R/C_6_R

  ENDIF

  RETURN
END SUBROUTINE INT_C1N6_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE INT_C1N6_RHS(imax,jkmax, ibc, dx, f,g)

  IMPLICIT NONE

  TINTEGER,                    INTENT(IN) :: imax, jkmax, ibc
  TREAL, DIMENSION(imax),      INTENT(IN) :: dx
  TREAL, DIMENSION(jkmax,imax),INTENT(IN) :: f
  TREAL, DIMENSION(jkmax,imax),INTENT(OUT):: g

! -------------------------------------------------------------------
  TINTEGER i
  TREAL c0102, c0103, c0104, c0106, c0203

! ###################################################################
  c0102 = C_1_R /C_2_R 
  c0103 = C_1_R /C_3_R 
  c0104 = C_1_R /C_4_R 
  c0106 = C_1_R /C_6_R 
  c0203 = C_2_R /C_3_R 

! -------------------------------------------------------------------
! Boundary conditions
! -------------------------------------------------------------------
  IF      ( ibc .EQ. 1 ) THEN ! array A22R
! BCs, see notes
     g(:,1) =-C_2_R*f(:,2)*dx(2) ! contribution to u'_1
!     g(:,2) = c0104*f(:,3)*dx(3) + c0102*f(:,2)*dx(2)
     g(:,2) = c0102*f(:,3)*dx(3) + c0203*f(:,2)*dx(2) ! fifth-order biased

! fourth-order centered
!     g(:,imax-1) = c0104*(f(:,imax-2)*dx(imax-2)+f(:,imax)*dx(imax)) + f(:,imax-1)*dx(imax-1)
! fifth-order biased
     g(:,imax-1) = c0102*f(:,imax-2)*dx(imax-2)+f(:,imax-1)*dx(imax-1)+c0106*f(:,imax)*dx(imax)
! third-order biased
     g(:,imax  ) = f(:,imax)*dx(imax  ) + C_2_R*f(:,imax-1)*dx(imax-1)

! -------------------------------------------------------------------
  ELSE IF ( ibc .EQ. 2 ) THEN ! array A11R
! third-order biased
     g(:,1) = f(:,1)*dx(1) + C_2_R*f(:,2)*dx(2)
! fourth-order centered
!     g(:,2) = c0104*(f(:,1)*dx(1)+f(:,3)*dx(3)) + f(:,2)*dx(2)
! fifth-order
     g(:,2) = c0106*f(:,1)*dx(1)+f(:,2)*dx(2)+c0102*f(:,3)*dx(3) ! fifth-order biased

! BCs, see notes
!     g(:,imax-1) = c0104*f(:,imax-2)*dx(imax-2) + c0102*f(:,imax-1)*dx(imax-1)
     g(:,imax-1) = c0102*f(:,imax-2)*dx(imax-2) + c0203*f(:,imax-1)*dx(imax-1)! fifth-order biased
     g(:,imax  ) =-C_2_R*f(:,imax-1)*dx(imax-1) ! contribution to u'_N
  ENDIF

! -------------------------------------------------------------------
! Interior points
! -------------------------------------------------------------------
! sixth-order
  DO i = 3,imax-2
     g(:,i) = f(:,i)*dx(i) + c0103*( f(:,i-1)*dx(i-1) + f(:,i+1)*dx(i+1) ) 
  ENDDO

  RETURN
END SUBROUTINE INT_C1N6_RHS
