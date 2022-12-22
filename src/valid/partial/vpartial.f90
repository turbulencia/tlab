!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/17 - J. Kostelecky
!#              Modified for pentadiagonal schemes
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate compact schemes.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

PROGRAM VPARTIAL

  USE TLAB_TYPES, ONLY : grid_dt
  USE TLAB_VARS,  ONLY : C1N6M_ALPHA
  use OPR_PARTIAL

  IMPLICIT NONE
 
#include "integers.h"
  
  TYPE(grid_dt)                      :: g
   
  TINTEGER                           :: jmax,kmax, i, l
  TINTEGER, PARAMETER                :: imax=128, len=1, inb_grid=57
   
  TREAL,    DIMENSION(imax,inb_grid) :: x
  TREAL,    DIMENSION(len,imax)      :: u, wrk3d
  TREAL,    DIMENSION(len,imax)      :: du1_a, du1_b, du1_c
  TREAL,    DIMENSION(len,imax)      :: du2_a
  !  TREAL,    DIMENSION(len,imax)      :: du2_n1, du2_n2, du2_n3
  TREAL,    DIMENSION(imax,7)        :: wrk1d
  TREAL,    DIMENSION(len)           :: wrk2d
  TREAL,    DIMENSION(len,2)         :: bcs
    TINTEGER bcs_aux(2,2)
  TREAL                              :: lambda, error, dummy
  TINTEGER                           :: test_type, ibc

! ###################################################################
! Initialize
  g%size     = imax 
  g%scale    = C_1_R
  g%uniform  = .TRUE.
  jmax       = 1
  kmax       = 1


! Valid stettings
  g%periodic = .FALSE.
  lambda     = 1 ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) lambda 
  test_type  = 2 
  ibc        = 3
  g%mode_fdm = FDM_COM6_JACOBIAN ! FDM_COM6_JACPENTA

  IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) C1N6M_ALPHA = 0.56
 
!  ###################################################################
   
  IF ( g%periodic ) THEN
    DO i = 1,imax
      x(i,1) = M_REAL(i-1)/M_REAL(imax)*g%scale
    ENDDO
  ELSE
    DO i = 1,imax
      x(i,1) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
    ENDDO
  ENDIF

  CALL FDM_INITIALIZE(x, g, wrk1d)
 
! Bcs
  bcs(1,1) = C_0_R
  bcs(1,2) = C_0_R
  bcs_aux = 0

! ###################################################################
! Define the function and analytic derivatives
  DO i = 1,imax
    DO l = 1,len
! single-mode
      u(l,i)     =                                 &
            SIN(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))!+C_PI_R/C_4_R)
      du1_a(l,i) = (C_2_R*C_PI_R/g%scale*lambda)    &
           *COS(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))!+C_PI_R/C_4_R)
      du2_a(l,i) =-(C_2_R*C_PI_R/g%scale*lambda)**2 &
           *u(l,i)
! Gaussian
      ! dummy = C_1_R / ( C_2_R*(g%scale/M_REAL(lambda*l))**2 )
      ! u(l,i)     = EXP(-dummy*(g%nodes(i)-x_0*g%scale)**2)
      ! du1_a(l,i) =-C_2_R *dummy *(g%nodes(i)-x_0*g%scale) *u(l,i)
      ! du2_a(l,i) =-C_2_R *dummy *(g%nodes(i)-x_0*g%scale) *du1_a(l,i) - C_2_R *dummy *u(l,i)
! Exponential
      ! u(l,i)     = EXP(g%nodes(i)/(g%scale/lambda))
      ! du1_a(l,i) = lambda/g%scale*u(l,i)
      ! du2_a(l,i) = lambda/g%scale*du1_a(l,i)
! delta-function
      ! u(i)     = MAX(C_0_R,C_2_R-M_REAL(i))
      ! du1_a(i) = C_0_R
      ! du2_a(i) = C_0_R
! hyperboic tangent
      ! u(l,i)     = lambda*LOG(C_1_R+EXP(g%nodes(i)/lambda))
      ! du1_a(l,i) = C_05_R*(C_1_R+TANH(C_05_R*g%nodes(i)/lambda))
      ! du2_a(l,i) = C_025_R/lambda/(COSH(C_05_R*g%nodes(i)/lambda))**2
! Polynomial
      ! dummy = C_4_R
      ! u(l,i)     =                       ( (g%scale-g%nodes(i)) /lambda)** dummy
      ! du1_a(l,i) = dummy                *( (g%scale-g%nodes(i)) /lambda)**(dummy-C_1_R)
      ! du2_a(l,i) = dummy *(dummy-C_1_R) *( (g%scale-g%nodes(i)) /lambda)**(dummy-C_2_R)
    ENDDO
  ENDDO

! ###################################################################
! Testing first-order derivatives
  
  IF ( test_type .EQ. 1 ) THEN 
    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_aux, g, u, du1_b, wrk3d, wrk2d,wrk3d)

! ! -------------------------------------------------------------------
! ! Testing second-order derivatives
! ! -------------------------------------------------------------------
!   ! Jacobian based
!     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g, u,     du2_n2, du1_n, wrk2d,wrk3d)
!     CALL OPR_PARTIAL_X(OPR_P1,    imax,jmax,kmax, bcs, g, du1_n, du2_n1, wrk3d, wrk2d,wrk3d)
!   ! Direct metrics
!     CALL FDM_C2N6ND_INITIALIZE(imax, x, wrk1d(1,1), wrk1d(1,4))
!     CALL TRIDFS(imax,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))  
!     CALL FDM_C2N6ND_RHS(imax,len, wrk1d(1,4), u, du2_n3)
!     CALL TRIDSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du2_n3)
   
! -------------------------------------------------------------------
  ELSEIF ( test_type .EQ. 2 ) THEN ! Testing new BCs routines
  
    IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
      CALL FDM_C1N6_BCS_LHS(imax,     ibc, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C1N6_BCS_RHS(imax,len, ibc,        u,du1_b)
    ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
      CALL FDM_C1N6M_BCS_LHS(imax,     ibc, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
      CALL FDM_C1N6M_BCS_RHS(imax,len, ibc,        u,du1_b)
    ENDIF

    IF      ( ibc .EQ. 0 ) THEN
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        CALL TRIDFS(imax,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du1_b)
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        CALL PENTADFS2(imax,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL PENTADSS2(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), du1_b)
      ENDIF

    ELSE IF ( ibc .EQ. 1 ) THEN
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        wrk1d(:,4) = C_0_R
        wrk1d(1,4) = C_1_R 
        wrk1d(2,4) = wrk1d(1,1) 
        wrk1d(3,4) = wrk1d(2,1)
        CALL TRIDFS(imax-1,     wrk1d(2,1),wrk1d(2,2),wrk1d(2,3))
        CALL TRIDSS(imax-1,len, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), du1_b(1,2))
        CALL TRIDSS(imax-1,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), wrk1d(2,4))
        bcs(:,1)   = du1_b(:,1)
        du1_b(:,1) = C_0_R
        DO i = 1,imax
          du1_b(:,i) = du1_b(:,i) + du1_a(:,1)*wrk1d(i,4) ! BCs
        ENDDO
        bcs(:,1) = bcs(:,1) + (wrk1d(1,2)*du1_b(:,1) + wrk1d(1,3)*du1_b(:,2))
        write(*,*) bcs(:,1), u(:,1)
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        wrk1d(:,6) = C_0_R
        wrk1d(1,6) = C_1_R 
        wrk1d(2,6) = wrk1d(1,2)
        wrk1d(3,6) = wrk1d(2,2)
        CALL PENTADFS2(imax-1,      wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))
        CALL PENTADSS2(imax-1, len, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), du1_b(1,2))
        CALL PENTADSS2(imax-1,  i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,6))
        bcs(:,1)   = du1_b(:,1)
        du1_b(:,1) = C_0_R
        DO i = 1,imax
          du1_b(:,i) = du1_b(:,i) + du1_a(:,1)*wrk1d(i,6) ! BCs
        ENDDO
        bcs(:,1) = bcs(:,1) + (wrk1d(1,3)*du1_b(:,1) + wrk1d(1,4)*du1_b(:,2))
        write(*,*) bcs(:,1), u(:,1)
      ENDIF

    ELSE IF ( ibc .EQ. 2 ) THEN
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        wrk1d(:     ,5) = C_0_R
        wrk1d(imax  ,5) = C_1_R
        wrk1d(imax-2,5) = wrk1d(imax-1,3)
        wrk1d(imax-1,5) = wrk1d(imax  ,3)
        CALL TRIDFS(imax-1,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(imax-1,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du1_b)
        CALL TRIDSS(imax-1,i1,  wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), wrk1d(1,5))
        bcs(:,2)      = du1_b(:,imax)
        du1_b(:,imax) = C_0_R
        DO i = 1,imax
           du1_b(:,i) = du1_b(:,i) + du1_a(:,imax)*wrk1d(i,5) ! BCs
        ENDDO
        bcs(:,2) = bcs(:,2) + (wrk1d(imax,1)*du1_b(:,imax-1) + wrk1d(imax,2)*du1_b(:,imax))
        write(*,*) bcs(:,2), u(:,imax)
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        wrk1d(:     ,6) = C_0_R
        wrk1d(imax  ,6) = C_1_R
        wrk1d(imax-2,6) = wrk1d(imax-1,4)
        wrk1d(imax-1,6) = wrk1d(imax  ,4)
        CALL PENTADFS2(imax-1,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL PENTADSS2(imax-1,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), du1_b)
        CALL PENTADSS2(imax-1,i1,  wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6))
        bcs(:,2)      = du1_b(:,imax)
        du1_b(:,imax) = C_0_R
        DO i = 1,imax
          du1_b(:,i) = du1_b(:,i) + du1_a(:,imax)*wrk1d(i,6) ! BCs
        ENDDO
        bcs(:,2) = bcs(:,2) + (wrk1d(imax,2)*du1_b(:,imax-1) + wrk1d(imax,3)*du1_b(:,imax))
        write(*,*) bcs(:,2), u(:,imax)
      ENDIF

    ELSE IF ( ibc .EQ. 3 ) THEN
      IF (g%mode_fdm .EQ. FDM_COM6_JACOBIAN) THEN
        wrk1d(:,4)      = C_0_R
        wrk1d(1,4)      = C_1_R 
        wrk1d(2,4)      = wrk1d(1,1)     
        wrk1d(3,4)      = wrk1d(2,1)
        wrk1d(:     ,5) = C_0_R
        wrk1d(imax  ,5) = C_1_R
        wrk1d(imax-2,5) = wrk1d(imax-1,3)
        wrk1d(imax-1,5) = wrk1d(imax  ,3)
        CALL TRIDFS(imax-2,     wrk1d(2,1),wrk1d(2,2),wrk1d(2,3))
        CALL TRIDSS(imax-2,len, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), du1_b(1,2))
        CALL TRIDSS(imax-2,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), wrk1d(2,4))
        CALL TRIDSS(imax-2,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), wrk1d(2,5))
        bcs(:,1)   = du1_b(:,1)
        du1_b(:,1) = C_0_R
        DO i = 1,imax
          du1_b(:,i) = du1_b(:,i) + du1_a(:,1)*wrk1d(i,4) ! BCs
        ENDDO
        bcs(:,1) = bcs(:,1) + (wrk1d(1,2)*du1_b(:,1) + wrk1d(1,3)*du1_b(:,2))
        write(*,*) bcs(:,1), u(:,1)
        bcs(:,2)      = du1_b(:,imax)
        du1_b(:,imax) = C_0_R
        DO i = 1,imax
          du1_b(:,i) = du1_b(:,i) + du1_a(:,imax)*wrk1d(i,5) ! BCs
        ENDDO
        bcs(:,2) = bcs(:,2) + (wrk1d(imax,1)*du1_b(:,imax-1) + wrk1d(imax,2)*du1_b(:,imax))
        write(*,*) bcs(:,2), u(:,imax)
      ELSEIF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) THEN
        wrk1d(:,6)      = C_0_R
        wrk1d(1,6)      = C_1_R 
        wrk1d(2,6)      = wrk1d(1,2)     
        wrk1d(3,6)      = wrk1d(2,2)
        wrk1d(:     ,7) = C_0_R
        wrk1d(imax  ,7) = C_1_R
        wrk1d(imax-2,7) = wrk1d(imax-1,4)
        wrk1d(imax-1,7) = wrk1d(imax  ,4)
        CALL PENTADFS2(imax-2,     wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))
        CALL PENTADSS2(imax-2,len, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), du1_b(1,2))
        CALL PENTADSS2(imax-2,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,6))
        CALL PENTADSS2(imax-2,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,7))
        bcs(:,1)   = du1_b(:,1)
        du1_b(:,1) = C_0_R
        DO i = 1,imax
          du1_b(:,i) = du1_b(:,i) + du1_a(:,1)*wrk1d(i,6) ! BCs
        ENDDO
        bcs(:,1) = bcs(:,1) + (wrk1d(1,3)*du1_b(:,1) + wrk1d(1,4)*du1_b(:,2))
        write(*,*) bcs(:,1), u(:,1)
        bcs(:,2)      = du1_b(:,imax)
        du1_b(:,imax) = C_0_R
        DO i = 1,imax
          du1_b(:,i) = du1_b(:,i) + du1_a(:,imax)*wrk1d(i,7) ! BCs
        ENDDO
        bcs(:,2) = bcs(:,2) + (wrk1d(imax,2)*du1_b(:,imax-1) + wrk1d(imax,3)*du1_b(:,imax))
        write(*,*) bcs(:,2), u(:,imax)
      ENDIF
    ENDIF
  ENDIF

! ###################################################################
! IO - Error and function values
  OPEN(20,file='partial.dat')
  error = C_0_R; dummy = C_0_R
  DO i = 1,imax
    DO l = 1,len
      ! Testing first-order derivatives
      WRITE(20,1000) g%nodes(i), u(l,i), du1_a(l,i), du1_b(l,i), du1_a(l,i) - du1_b(l,i)
      du1_c(l,i) = ABS(du1_a(l,i) - du1_b(l,i))
      dummy = dummy + du1_a(l,i) * du1_a(l,i)
      error = error + du1_c(l,i) * du1_c(l,i)
      ! ! Testing second-order derivatives 
      ! WRITE(20,1000) g%nodes(i), u(l,i), du2_a(l,i), du2_n2(l,i), du2_a(l,i)-du2_n2(l,i)
      ! du1_c(l,i)= ABS(du2_a(l,i)-du2_n2(l,i))
      ! dummy = dummy + du2_a(l,i)*du2_a(l,i)
      ! error = error + du1_c(l,i) * du1_c(l,i)
    ENDDO
  ENDDO
  CLOSE(20)

  WRITE(*,*) 'Solution L2-norm ...........:', SQRT(g%jac(1,1)*dummy) / M_REAL(len)
  IF ( dummy .EQ. C_0_R ) STOP
  WRITE(*,*) 'Relative Error L2-norm .....:', SQRT(g%jac(1,1)*error)  /MAXVAL(ABS(du1_a))
  WRITE(*,*) 'Relative Error Linf-norm ...:', MAXVAL(du1_c(1,1:imax)) /MAXVAL(ABS(du1_a))

  STOP

1000 FORMAT(5(1x,e12.5))

END PROGRAM VPARTIAL