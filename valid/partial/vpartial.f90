#include "types.h"
#include "dns_const.h"

PROGRAM VPARTIAL

  USE DNS_TYPES, ONLY : grid_dt
  IMPLICIT NONE

#include "integers.h"

  TYPE(grid_dt) :: g
  TINTEGER imax,jmax,kmax, i, l, ibc, len, inb_grid
  PARAMETER(imax=128, len=1, inb_grid=3+4*3+4*3+1*5)
  TREAL lambda, x_0
  TREAL, DIMENSION(imax,inb_grid) :: x
  TREAL, DIMENSION(len,imax)      :: u, wrk3d
  TREAL, DIMENSION(len,imax)      :: du1_a, du1_n
  TREAL, DIMENSION(len,imax)      :: du2_a, du2_n1, du2_n2, du2_n3
  TREAL wrk1d(imax,7), wrk2d(len), bcs(len,2)
  TREAL error, dummy

  TINTEGER type

! ###################################################################
  g%size     = imax
  g%scale    = C_1_R
  g%mode_fdm = FDM_COM6_JACOBIAN
  g%uniform  = .TRUE.
  jmax = 1
  kmax = 1

!  WRITE(*,*) 'Periodic (0) or nonperiodic (1) case ?'
  g%periodic = .FALSE.

  WRITE(*,*) 'Parameter ?'
  READ(*,*) lambda

  IF ( g%periodic ) THEN
     DO i = 1,imax
        x(i,1) = M_REAL(i-1)/M_REAL(imax)*g%scale
     ENDDO

  ELSE
     DO i = 1,imax
        x(i,1) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
     ENDDO
  ENDIF
  ! OPEN(21,file='grid.dat')
  ! DO i = 1,imax
  !    READ(21,*) x(i,1)
  ! ENDDO
  ! g%scale=x(imax,1)-x(1,1)
  ! CLOSE(21)
  ! g%uniform  = .FALSE.

  CALL FDM_INITIALIZE(x, g, wrk1d)

  x_0 = C_05_R

! ###################################################################
! Define the function
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
!  u(i)     = MAX(C_0_R,C_2_R-M_REAL(i))
!  du1_a(i) = C_0_R
!  du2_a(i) = C_0_R
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
  type = 1

! ###################################################################
  IF ( type .EQ. 1 ) THEN ! Testing second-order derivatives
! Jacobian based
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g, u,     du2_n2, du1_n, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1,    imax,jmax,kmax, bcs, g, du1_n, du2_n1, wrk3d, wrk2d,wrk3d)

! Direct metrics
     CALL FDM_C2N6ND_INITIALIZE(imax, x, wrk1d(1,1), wrk1d(1,4))
     CALL TRIDFS(imax,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))

     CALL FDM_C2N6ND_RHS(imax,len, wrk1d(1,4), u, du2_n3)
     CALL TRIDSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du2_n3)

  ENDIF

! -------------------------------------------------------------------
  IF ( type .EQ. 2 ) THEN ! Testing new BCs routines
     ibc = 0
     CALL FDM_C1N6_BCS_LHS(imax,     ibc, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
     CALL FDM_C1N6_BCS_RHS(imax,len, ibc,        u,du1_n)

     IF      ( ibc .EQ. 0 ) THEN
        CALL TRIDFS(imax,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du1_n)

     ELSE IF ( ibc .EQ. 1 ) THEN
        wrk1d(:,4) = C_0_R; wrk1d(1,4) =  C_1_R; wrk1d(2,4) =  wrk1d(1,1); wrk1d(3,4) =  wrk1d(2,1)
        CALL TRIDFS(imax-1,     wrk1d(2,1),wrk1d(2,2),wrk1d(2,3))
        CALL TRIDSS(imax-1,len, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), du1_n(1,2))
        CALL TRIDSS(imax-1,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), wrk1d(2,4))
        bcs(:,1) = du1_n(:,1); du1_n(:,1) = C_0_R
        DO i = 1,imax
           du1_n(:,i) = du1_n(:,i) + du1_a(:,1)*wrk1d(i,4) ! BCs
        ENDDO
        bcs(:,1) = bcs(:,1) + (wrk1d(1,2)*du1_n(:,1) + wrk1d(1,3)*du1_n(:,2))
        print*, bcs(:,1), u(:,1)

     ELSE IF ( ibc .EQ. 2 ) THEN
        wrk1d(:,5) = C_0_R; wrk1d(imax,5) = C_1_R; wrk1d(imax-2,5) = wrk1d(imax-1,3); wrk1d(imax-1,5) = wrk1d(imax,3)
        CALL TRIDFS(imax-1,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(imax-1,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du1_n)
        CALL TRIDSS(imax-1,i1,  wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), wrk1d(1,5))
        bcs(:,2) = du1_n(:,imax); du1_n(:,imax) = C_0_R
        DO i = 1,imax
           du1_n(:,i) = du1_n(:,i) + du1_a(:,imax)*wrk1d(i,5) ! BCs
        ENDDO
        bcs(:,2) = bcs(:,2) + (wrk1d(imax,1)*du1_n(:,imax-1) + wrk1d(imax,2)*du1_n(:,imax))
        print*, bcs(:,2), u(:,imax)

     ELSE IF ( ibc .EQ. 3 ) THEN
        wrk1d(:,4) = C_0_R; wrk1d(1,4)    = C_1_R; wrk1d(2,4)      = wrk1d(1,1);      wrk1d(3,4)      = wrk1d(2,1)
        wrk1d(:,5) = C_0_R; wrk1d(imax,5) = C_1_R; wrk1d(imax-2,5) = wrk1d(imax-1,3); wrk1d(imax-1,5) = wrk1d(imax,3)
        CALL TRIDFS(imax-2,     wrk1d(2,1),wrk1d(2,2),wrk1d(2,3))
        CALL TRIDSS(imax-2,len, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), du1_n(1,2))
        CALL TRIDSS(imax-2,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), wrk1d(2,4))
        CALL TRIDSS(imax-2,i1,  wrk1d(2,1),wrk1d(2,2),wrk1d(2,3), wrk1d(2,5))

        bcs(:,1) = du1_n(:,1); du1_n(:,1) = C_0_R
        DO i = 1,imax
           du1_n(:,i) = du1_n(:,i) + du1_a(:,1)*wrk1d(i,4) ! BCs
        ENDDO
        bcs(:,1) = bcs(:,1) + (wrk1d(1,2)*du1_n(:,1) + wrk1d(1,3)*du1_n(:,2))
        print*, bcs(:,1), u(:,1)

        bcs(:,2) = du1_n(:,imax); du1_n(:,imax) = C_0_R
        DO i = 1,imax
           du1_n(:,i) = du1_n(:,i) + du1_a(:,imax)*wrk1d(i,5) ! BCs
        ENDDO
        bcs(:,2) = bcs(:,2) + (wrk1d(imax,1)*du1_n(:,imax-1) + wrk1d(imax,2)*du1_n(:,imax))
        print*, bcs(:,2), u(:,imax)

     ENDIF

  ENDIF

! ###################################################################
  OPEN(20,file='partial.dat')
  error = C_0_R; dummy = C_0_R
  DO i = 1,imax
     DO l = 1,len
        ! WRITE(20,1000) g%nodes(i), u(l,i), du1_a(l,i), du1_n(l,i), du1_a(l,i)-du1_n(l,i)
        ! du1_n(l,i)= ABS(du1_a(l,i)-du1_n(l,i))
        ! dummy   = dummy   + du1_a(l,i)*du1_a(l,i)
        WRITE(20,1000) g%nodes(i), u(l,i), du2_a(l,i), du2_n2(l,i), du2_a(l,i)-du2_n2(l,i)
        du1_n(l,i)= ABS(du2_a(l,i)-du2_n2(l,i))
        dummy = dummy + du2_a(l,i)*du2_a(l,i)

        error = error + du1_n(l,i)*du1_n(l,i)
     ENDDO
  ENDDO
  CLOSE(20)

  WRITE(*,*) 'Solution L2-norm ...........:', SQRT(g%jac(1,1)*dummy) / M_REAL(len)
  IF ( dummy .EQ. C_0_R ) STOP

  WRITE(*,*) 'Relative Error L2-norm .....:', SQRT(g%jac(1,1)*error)  /MAXVAL(ABS(du2_a))
  WRITE(*,*) 'Relative Error Linf-norm ...:', MAXVAL(du1_n(1,1:imax)) /MAXVAL(ABS(du2_a))

  STOP

1000 FORMAT(5(1x,e12.5))

END PROGRAM VPARTIAL
