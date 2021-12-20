!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2010/12/20 - J. Kostelecky
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
!# The following IVP is solved
!#
!#     u'_i + \lamba u_i = h_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!# The system of N-1 eqns:
!# 
!#                    (B + \lambda A)u = Ah = l
!#
!# is established in this routine, giving diagonals a-g and l (see notes).
!#
!# Carpenter, Gottlieb and Aberbanel, JCP, 1993, study the effect of 
!# boundary points on stability. Scheme 3-5-2xtri6c--penta6c--2xtri6c-5-3 
!# is implmented.
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
!# h           Out   forcing term for the exponential
!#
!########################################################################
#include "types.h"

!########################################################################
!Left-hand side; heptadiagonal matrix of the linear system and h
!########################################################################
SUBROUTINE INT_C1N6M_LHS_E(imax, ibc, dx, lambda, a,b,c,d,e,f,g, h)

  USE TLAB_VARS, ONLY : C1N6M_ALPHA, C1N6M_BETA
  USE TLAB_VARS, ONLY : C1N6M_AD2, C1N6M_BD4, C1N6M_CD6

  IMPLICIT NONE

  TREAL lambda
  TINTEGER,               INTENT(IN)  :: imax, ibc
  TREAL, DIMENSION(imax), INTENT(IN)  :: dx
  TREAL, DIMENSION(imax), INTENT(OUT) :: a,b,c,d,e,f,g,h

! -------------------------------------------------------------------
  TINTEGER                            :: i
  TREAL                               :: c0136, c1418, c0103, c0104

! ###################################################################
  c0136 = C_1_R /C_36_R
  c1418 = C_14_R/C_18_R
  c0103 = C_1_R /C_3_R 
  c0104 = C_1_R /C_4_R 

! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
  a(1     ) =   C_0_R 
  b(1     ) =   C_0_R 
  c(1     ) =   C_0_R 
  d(1     ) = - C_5_R/C_2_R + lambda       *dx(1)
  e(1     ) =   C_2_R       + lambda*C_2_R *dx(2)
  f(1     ) =   C_05_R
  g(1     ) =   C_0_R
! fifth-order biased
  a(2     ) =   C_0_R 
  b(2     ) =   C_0_R 
  c(2     ) = - C_5_R/C_9_R + lambda*C_1_R/C_6_R *dx(1)
  d(2     ) = - C_1_R/C_2_R + lambda             *dx(2)
  e(2     ) =   C_1_R       + lambda*C_1_R/C_2_R *dx(3)
  f(2     ) =   C_1_R/C_18_R 
  g(2     ) =   C_0_R 
! sixth-order centered (alpha=1/3)
  DO i = 3,4 ! DO i = 3,3
    a(i   ) =   C_0_R
    b(i   ) = - c0136
    c(i   ) = - c1418       + lambda*c0103  *dx(i-1)
    d(i   ) =   C_0_R       + lambda        *dx(i  )
    e(i   ) =   c1418       + lambda*c0103  *dx(i+1)
    f(i   ) =   c0136
    g(i   ) =   C_0_R
  ENDDO
! sixth-order modified centered
  DO i = 5,imax-4 ! DO i = 4,imax-3
    a(i   ) = - C1N6M_CD6
    b(i   ) = - C1N6M_BD4   + lambda* C1N6M_BETA  *dx(i-2)
    c(i   ) = - C1N6M_AD2   + lambda* C1N6M_ALPHA *dx(i-1)
    d(i   ) =   C_0_R       + lambda              *dx(i  )
    e(i   ) =   C1N6M_AD2   + lambda* C1N6M_ALPHA *dx(i+1)
    f(i   ) =   C1N6M_BD4   + lambda* C1N6M_BETA  *dx(i+2)
    g(i   ) =   C1N6M_CD6
  ENDDO
! sixth-order centered (alpha=1/3)
  DO i = imax-3,imax-2 ! DO i = imax-2,imax-2
    a(i   ) =   C_0_R
    b(i   ) = - c0136
    c(i   ) = - c1418       + lambda*c0103  *dx(i-1)
    d(i   ) =   C_0_R       + lambda        *dx(i  )
    e(i   ) =   c1418       + lambda*c0103  *dx(i+1)
    f(i   ) =   c0136
    g(i   ) =   C_0_R
  ENDDO
! fifth-order biased
  a(imax-1) =   C_0_R
  b(imax-1) = - C_1_R/C_18_R
  c(imax-1) = - C_1_R       + lambda*C_1_R/C_2_R *dx(imax-2)
  d(imax-1) =   C_1_R/C_2_R + lambda             *dx(imax-1)
  e(imax-1) =   C_5_R/C_9_R + lambda*C_1_R/C_6_R *dx(imax  )
  f(imax-1) =   C_0_R 
  g(imax-1) =   C_0_R 
! third-order biased
  a(imax  ) =   C_0_R      
  b(imax  ) = - C_05_R      
  c(imax  ) = - C_2_R       + lambda*C_2_R *dx(imax-1)
  d(imax  ) =   C_5_R/C_2_R + lambda       *dx(imax  )
  e(imax  ) =   C_0_R 
  f(imax  ) =   C_0_R 
  g(imax  ) =   C_0_R 
    
! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
  IF      ( ibc .EQ. 1 ) THEN
    d(1     ) = - C_5_R/C_2_R                                   ! array B22R
    d(2     ) = - C_5_R /C_6_R  + lambda*C_2_R/C_3_R*dx(2)      ! B22R + lambda A22R fifth-order biased
    e(2     ) =   C_11_R/C_12_R + lambda*C_1_R/C_2_R*dx(3)     
    f(imax  ) =   C_1_R/C_6_R                                   ! element A(imax-1,imax)

  ELSE IF ( ibc .EQ. 2 ) THEN
    b(1     ) =   C_1_R/C_6_R                                   ! element A(2,1)
    c(imax-1) = - C_11_R/C_12_R + lambda*C_1_R/C_2_R*dx(imax-2) ! B11R + lambda A11R fifth-order biased
    d(imax-1) =   C_5_R /C_6_R  + lambda*C_2_R/C_3_R*dx(imax-1)
    d(imax  ) =   C_5_R/C_2_R                                   ! array B11R

  ENDIF

! -------------------------------------------------------------------
! Setting the RHS for the null space (with minus sign)
! -------------------------------------------------------------------
  h = C_0_R
  IF      ( ibc .EQ. 1 ) THEN
    h(1     ) =   C_1_R       ! normalization
    h(2     ) =   c0136*C_5_R ! fifth-order biased
    h(3     ) =   c0136
  ELSE IF ( ibc .EQ. 2 ) THEN
    h(imax-2) = - c0136       ! fifth-order biased
    h(imax-1) = - c0136*C_5_R
    h(imax  ) =   C_1_R       ! normalization
  ENDIF

  RETURN
END SUBROUTINE INT_C1N6M_LHS_E

!########################################################################
!Left-hand side case \lambda=0; pentadiagonal matrix of the linear system
!########################################################################
SUBROUTINE INT_C1N6M_LHS(imax, ibc, a,b,c,d,e,f,g)

  USE TLAB_VARS, ONLY : C1N6M_AD2, C1N6M_BD4, C1N6M_CD6

  IMPLICIT NONE

  TINTEGER,               INTENT(IN) :: imax, ibc
  TREAL, DIMENSION(imax), INTENT(OUT):: a,b,c,d,e,f,g

! -------------------------------------------------------------------
  TREAL                              :: c0136, c1418
  TINTEGER                           :: i

! ###################################################################
  c0136 = C_1_R /C_36_R
  c1418 = C_14_R/C_18_R

! -------------------------------------------------------------------
! Define diagonals of heptadiagonal system
! -------------------------------------------------------------------
! third-order biased
  a(1       ) =   C_0_R 
  b(1       ) =   C_0_R 
  c(1       ) =   C_0_R 
  d(1       ) = - C_5_R/C_2_R
  e(1       ) =   C_2_R
  f(1       ) =   C_05_R      
  g(1       ) =   C_0_R      
! fifth-order biased
  a(2       ) =   C_0_R 
  b(2       ) =   C_0_R 
  c(2       ) = - C_5_R/C_9_R
  d(2       ) = - C_1_R/C_2_R
  e(2       ) =   C_1_R
  f(2       ) =   C_1_R/C_18_R
  g(2       ) =   C_0_R
! sixth-order centered (alpha=1/3)
  DO i = 3,4 ! DO i = 3,3
    a(i     ) =   C_0_R
    b(i     ) = - c0136
    c(i     ) = - c1418
    d(i     ) =   C_0_R
    e(i     ) =   c1418
    f(i     ) =   c0136
    g(i     ) =   C_0_R
  ENDDO
! sixth-order modified centered
  DO i = 5,imax-4 ! DO i = 4,imax-3
    a(i     ) = - C1N6M_CD6
    b(i     ) = - C1N6M_BD4
    c(i     ) = - C1N6M_AD2
    d(i     ) =   C_0_R
    e(i     ) =   C1N6M_AD2
    f(i     ) =   C1N6M_BD4
    g(i     ) =   C1N6M_CD6
  ENDDO
! sixth-order centered (alpha=1/3)
  DO i = imax-3,imax-2 ! DO i = imax-2,imax-2
    a(i     ) =   C_0_R
    b(i     ) = - c0136
    c(i     ) = - c1418
    d(i     ) =   C_0_R
    e(i     ) =   c1418
    f(i     ) =   c0136
    g(i     ) =   C_0_R
  ENDDO
! fifth-order biased
  a(  imax-1) =   C_0_R
  b(  imax-1) = - C_1_R/C_18_R
  c(  imax-1) = - C_1_R
  d(  imax-1) =   C_1_R/C_2_R      
  e(  imax-1) =   C_5_R/C_9_R
  f(  imax-1) =   C_0_R 
  g(  imax-1) =   C_0_R 
! third-order biased
  a(  imax  ) =   C_0_R      
  b(  imax  ) = - C_05_R      
  c(  imax  ) = - C_2_R       
  d(  imax  ) =   C_5_R/C_2_R 
  e(  imax  ) =   C_0_R 
  f(  imax  ) =   C_0_R 
  g(  imax  ) =   C_0_R 
     
! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
  IF      ( ibc .EQ. 1 ) THEN ! array B22R
    d(2     ) = - C_5_R /C_6_R  ! fifth-order biased
    e(2     ) =   C_11_R/C_12_R
    f(imax  ) =   C_1_R/C_6_R
         
  ELSE IF ( ibc .EQ. 2 ) THEN ! array B11R
    c(imax-1) = - C_11_R/C_12_R ! fifth-order biased
    d(imax-1) =   C_5_R /C_6_R
    b(1     ) =   C_1_R/C_6_R
    
  ENDIF

  RETURN
END SUBROUTINE INT_C1N6M_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE INT_C1N6M_RHS(imax,jkmax, ibc, dx, h,l)

  USE TLAB_VARS, ONLY : C1N6M_ALPHA, C1N6M_BETA

  IMPLICIT NONE

  TINTEGER,                    INTENT(IN) :: imax, jkmax, ibc
  TREAL, DIMENSION(imax),      INTENT(IN) :: dx
  TREAL, DIMENSION(jkmax,imax),INTENT(IN) :: h
  TREAL, DIMENSION(jkmax,imax),INTENT(OUT):: l

! -------------------------------------------------------------------
  TINTEGER                                :: i
  TREAL                                   :: c0102, c0103, c0104, c0106, c0203

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
    l(:,1     ) = - C_2_R*h(:,2)*dx(2)                      ! contribution to u'_1
    l(:,2     ) =   c0102*h(:,3)*dx(3) + c0203*h(:,2)*dx(2) ! fifth-order biased
! fifth-order biased
    l(:,imax-1) =   c0102*h(:,imax-2)*dx(imax-2)+h(:,imax-1)*dx(imax-1)+c0106*h(:,imax)*dx(imax)
! third-order biased
    l(:,imax  ) =   h(:,imax)*dx(imax  ) + C_2_R*h(:,imax-1)*dx(imax-1)

! -------------------------------------------------------------------
  ELSE IF ( ibc .EQ. 2 ) THEN ! array A11R
! third-order biased
    l(:,1     ) =   h(:,1)*dx(1) + C_2_R*h(:,2)*dx(2)
! fifth-order
    l(:,2     ) =   c0106*h(:,1)*dx(1)+h(:,2)*dx(2)+c0102*h(:,3)*dx(3)          ! fifth-order biased
! BCs, see notes
    l(:,imax-1) =   c0102*h(:,imax-2)*dx(imax-2) + c0203*h(:,imax-1)*dx(imax-1) ! fifth-order biased
    l(:,imax  ) = - C_2_R*h(:,imax-1)*dx(imax-1)                                ! contribution to u'_N
  ENDIF

! -------------------------------------------------------------------
! Interior points
! -------------------------------------------------------------------
! 6th-order centered with alpha=(1/3)
  DO i = 3,4 ! DO i = 3,3
  l(:,i       ) =   h(:,i)*dx(i) + c0103*( h(:,i-1)*dx(i-1)  + h(:,i+1)*dx(i+1) ) 
  ENDDO
  !
  DO i = imax-3,imax-2 ! DO i = imax-2,imax-2
  l(:,i       ) =   h(:,i)*dx(i) + c0103*( h(:,i-1)*dx(i-1)  + h(:,i+1)*dx(i+1) ) 
  ENDDO

! sixth-order modified centered
  DO i = 5,imax-4 ! DO i = 4,imax-3
    l(:,i     ) =                   h(:,i  )*dx(i  )                       + &
                    C1N6M_ALPHA * ( h(:,i-1)*dx(i-1)  + h(:,i+1)*dx(i+1) ) + &
                    C1N6M_BETA  * ( h(:,i-2)*dx(i-2)  + h(:,i+2)*dx(i+2) )
  ENDDO

  RETURN
END SUBROUTINE INT_C1N6M_RHS