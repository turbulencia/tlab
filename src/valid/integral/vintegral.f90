!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/20 - J. Kostelecky
!#              Modified for pentadiagonal schemes
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

PROGRAM VINTEGRAL

  USE TLAB_TYPES, ONLY : grid_dt
  USE TLAB_VARS,  ONLY : C1N6M_ALPHA
  use OPR_PARTIAL

  IMPLICIT NONE

#include "integers.h"

  TYPE(grid_dt)                      :: g
  TINTEGER                           :: jmax,kmax, i
  TINTEGER, PARAMETER                :: imax=128, inb_grid=57
  TREAL,    DIMENSION(imax,inb_grid) :: x
  TREAL,    DIMENSION(imax         ) :: wrk2d, wrk3d
  TREAL,    DIMENSION(imax         ) :: u, w_n, f
  TREAL,    DIMENSION(imax         ) :: du1_a, dw1_n, du2_a
  TREAL,    DIMENSION(imax,8       ) :: wrk1d
  integer,    DIMENSION(2   ,2       ) :: bcs
  TREAL                              :: lambda, error, sol, dummy, wk, x_0
  TINTEGER                           :: test_type, ibc
  TINTEGER                           :: imin_loc, imax_loc

! ###################################################################
! Initialize
  g%size     = imax 
  g%scale    = C_1_R
  g%uniform  = .TRUE.
  jmax       = 1
  kmax       = 1
  bcs        = 0

! Valid stettings
  g%periodic = .TRUE.
  wk         = 1 ! WRITE(*,*) 'Wavenumber ?'; READ(*,*) wk
  lambda     = 1 ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) lambda 
  test_type  = 1
  g%mode_fdm = FDM_COM6_JACOBIAN ! FDM_COM6_JACPENTA

  IF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) C1N6M_ALPHA = 0.56

! ###################################################################

  DO i = 1,imax
     x(i,1) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
  ENDDO

  CALL FDM_INITIALIZE(x, g, wrk1d)
  
  x_0 = C_075_R
  
  wrk1d = C_0_R; wrk2d = C_0_R; wrk3d = C_0_R

! ###################################################################
! Define the function f and analytic derivatives
  DO i = 1,imax
! single-mode
    u(i)     = SIN(C_2_R*C_PI_R/g%scale*wk*g%nodes(i)+C_PI_R/C_4_R)
    du1_a(i) = (C_2_R*C_PI_R/g%scale*wk)&
              *COS(C_2_R*C_PI_R/g%scale*wk*g%nodes(i)+C_PI_R/C_4_R)
    ! u(i)     =              COS(C_2_R*C_PI_R*g%nodes(i))
    ! du1_a(i) =-C_PI_R*C_2_R*SIN(C_2_R*C_PI_R*g%nodes(i))
! Gaussian
    ! u(i)     = EXP(-(g%nodes(i)-x_0*g%scale)**2/(C_2_R*(g%scale/wk)**2))
    ! du1_a(i) =-(g%nodes(i)-x_0*g%scale)/(g%scale/wk)**2*u(i)
    ! v(i)     = EXP(-(g%nodes(i)-x_0*C_05_R*g%scale)**2/(C_2_R*(g%scale/wk)**2))
    ! dv1_a(i) =-(g%nodes(i)-x_0*C_05_R*g%scale)/(g%scale/wk)**2*v(i)
! exponential
    ! u(i) = EXP(-g%nodes(i)*lambda)
    ! du1_a(i) = -lambda*u(i)
    ! du2_a(i) =  lambda*lambda*u(i)
    ! v(i) = EXP(g%nodes(i)*x_0/g%scale)
    ! dv1_a(i) = x_0/g%scale*v(i)
! step
    ! u(i) = MAX(C_0_R,(g%nodes(i)-g%nodes(imax/2))*x_0)
    ! du1_a(i) = (C_1_R+SIGN(C_1_R,g%nodes(i)-g%nodes(imax/2)))*C_05_R*x_0
! tanh
    ! u(i) = x_0 * LOG( C_1_R + EXP( (g%nodes(i)-g%nodes(imax/2))/x_0 ) )
    ! du1_a(i) = C_05_R*( C_1_R + TANH( C_05_R*(g%nodes(i)-g%nodes(imax/2))/x_0 ) )
! polynomial
    ! u(i)     =         g%nodes(i)** lambda
    ! du1_a(i) = lambda*(g%nodes(i)**(lambda-C_1_R))
! zero
    ! u(i) = C_0_R
    ! du1_a(i) = C_0_R
  ENDDO

! ###################################################################
! Integral
! ###################################################################
  IF ( test_type .EQ. 0 ) THEN

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g, u, f, wrk3d, wrk2d,wrk3d) ! f = du1_a

    ibc = 2

    IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
      CALL INT_C1N6_LHS(imax,    ibc,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
      CALL INT_C1N6_RHS(imax,i1, ibc, g%jac, f,w_n)
    ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
      CALL INT_C1N6M_LHS(imax,    ibc,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7))
      CALL INT_C1N6M_RHS(imax,i1, ibc, g%jac, f,w_n)
    ENDIF 

    IF      ( ibc .EQ. 1 ) THEN ! at the bottom
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        CALL PENTADFS(imax-1,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))
        CALL PENTADSS(imax-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), w_n(2))
        w_n(1) = C_0_R; w_n   = w_n + u(1)    ! BCs
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        CALL HEPTADFS(imax-1,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),wrk1d(2,6),wrk1d(2,7))
        CALL HEPTADSS(imax-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),wrk1d(2,6),wrk1d(2,7), w_n(2))
        w_n(1) = C_0_R; w_n   = w_n + u(1)    ! BCs
      ENDIF 

    ELSE IF ( ibc .EQ. 2 ) THEN ! at the top
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        CALL PENTADFS(imax-1,    wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL PENTADSS(imax-1,i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), w_n(1))
        w_n(imax) = C_0_R; w_n   = w_n + u(imax) ! BCs
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        CALL HEPTADFS(imax-1,    wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7))
        CALL HEPTADSS(imax-1,i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7), w_n(1))
        w_n(imax) = C_0_R; w_n   = w_n + u(imax) ! BCs
      ENDIF 
    ENDIF

! ###################################################################
! First order equation
! ###################################################################
  ELSE IF ( test_type .EQ. 1 ) THEN

    f     = du1_a + lambda*u
    
    ibc   = 2

    IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
      CALL INT_C1N6_LHS_E( imax,    ibc, g%jac, lambda, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6))
      CALL INT_C1N6_RHS  ( imax,i1, ibc, g%jac,         f,w_n)
    ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
      CALL INT_C1N6M_LHS_E(imax,    ibc, g%jac, lambda, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7), wrk1d(1,8))
      CALL INT_C1N6M_RHS  (imax,i1, ibc, g%jac,         f,w_n)
    ENDIF 

    IF      ( ibc .EQ. 1 ) THEN
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        CALL PENTADFS(imax-1,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))
        CALL PENTADSS(imax-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), w_n(2))
        CALL PENTADSS(imax-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,6))
        dummy = w_n(1); w_n(1) = C_0_R
        w_n   = w_n + u(1)*wrk1d(1:imax,6) ! BCs
        dummy =(dummy+ wrk1d(1,3)*w_n(1)+ wrk1d(1,4)*w_n(2)+ wrk1d(1,5)*w_n(3))/g%jac(1,1)
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        CALL HEPTADFS(imax-1,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),wrk1d(2,6),wrk1d(2,7))
        CALL HEPTADSS(imax-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),wrk1d(2,6),wrk1d(2,7), w_n(2))
        CALL HEPTADSS(imax-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),wrk1d(2,6),wrk1d(2,7), wrk1d(2,8))
        dummy = w_n(1); w_n(1) = C_0_R
        w_n   = w_n + u(1)*wrk1d(1:imax,8) ! BCs
        dummy =(dummy+ wrk1d(1,4)*w_n(1)+ wrk1d(1,5)*w_n(2)+ wrk1d(1,6)*w_n(3))/g%jac(1,1)
      ENDIF 

    ELSE IF ( ibc .EQ. 2 ) THEN
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        CALL PENTADFS(imax-1,    wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL PENTADSS(imax-1,i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), w_n(1))
        CALL PENTADSS(imax-1,i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6))
        dummy = w_n(imax); w_n(imax) = C_0_R
        w_n   = w_n + u(imax)*wrk1d(1:imax,6) ! BCs
        dummy =(dummy+ wrk1d(imax,1)*w_n(imax-2)+ wrk1d(imax,2)*w_n(imax-1)+ wrk1d(imax,3)*w_n(imax))/g%jac(imax,1)
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        CALL HEPTADFS(imax-1,    wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7))
        CALL HEPTADSS(imax-1,i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7), w_n(1))
        CALL HEPTADSS(imax-1,i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7), wrk1d(1,8))
        dummy = w_n(imax); w_n(imax) = C_0_R
        w_n   = w_n + u(imax)*wrk1d(1:imax,8) ! BCs
        dummy =(dummy+ wrk1d(imax,2)*w_n(imax-2)+ wrk1d(imax,3)*w_n(imax-1)+ wrk1d(imax,4)*w_n(imax))/g%jac(imax,1)
      ENDIF 
    ENDIF

    IF ( ibc .EQ. 1 ) THEN
      WRITE(*,*) dummy, dw1_n(1   )
    ELSE
      WRITE(*,*) dummy, dw1_n(imax)
    ENDIF
    
    f = u

! ###################################################################
! Second order equation
! ###################################################################
  ELSE IF ( test_type .EQ. 2 ) THEN

    CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g, u, f, wrk3d, wrk2d,wrk3d)
     
    dummy = lambda*lambda
    f     = du2_a - dummy*u

    CALL INT_C2N6_LHS_E(imax,    g%jac, dummy, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6),wrk1d(1,7))
    CALL PENTADFS(imax-2,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))
    CALL PENTADSS(imax-2,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,6))
    CALL PENTADSS(imax-2,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,7))
    CALL INT_C2N6_RHS  (imax,i1, g%jac,         f,w_n)
    CALL PENTADSS(imax-2,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), w_n(2))
    w_n(:) = w_n(:) + u(1)*wrk1d(:,6) + u(imax)*wrk1d(:,7)
  ENDIF

! ###################################################################
! IO - Error and function values
  OPEN(20,file='integral.dat')
  error = C_0_R
  sol   = C_0_R
  imin_loc = 1; imax_loc = imax
  DO i = imin_loc,imax_loc
    WRITE(20,1000) g%nodes(i), u(i), w_n(i), u(i)-w_n(i)
    w_n(i)= ABS(u(i)-w_n(i))
    error = error + w_n(i)*w_n(i)
    sol   = sol   + u(i)*u(i)
  ENDDO
  CLOSE(20)

  WRITE(*,*) 'Solution L2-norm ......:', sqrt(g%jac(1,1)*sol)
  WRITE(*,*) 'Error L2-norm .........:', sqrt(g%jac(1,1)*error)
  WRITE(*,*) 'Error Linf-norm .......:', MAXVAL(w_n(1:imax))
  WRITE(*,*) 'Relative error ........:', sqrt(error)/sqrt(sol)
  WRITE(*,*) 'Derivative overshoot ..:', MINVAL(dw1_n(1:imax))

  STOP

1000 FORMAT(6(1x,e17.10e3))

END PROGRAM VINTEGRAL