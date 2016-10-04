#include "types.h"

PROGRAM VPARTIAL
  
  USE DNS_GLOBAL, ONLY : imax_total, inb_scal

  IMPLICIT NONE
  
#include "integers.h"
  
  TINTEGER imode_fdm, imax,jmax,kmax, i, l, i1bc, iunif, ibc, len, inb_grid
  PARAMETER(imax=3200, len=1, inb_grid=2+4*3+4*3+1*5)
  TREAL scalex, lambda, x_0
  TREAL, DIMENSION(imax)          :: x
  TREAL, DIMENSION(imax,inb_grid) :: dx
  TREAL, DIMENSION(len,imax)      :: u, wrk3d
  TREAL, DIMENSION(len,imax)      :: du1_a, du1_n
  TREAL, DIMENSION(len,imax)      :: du2_a, du2_n1, du2_n2, du2_n3
  TREAL wrk1d(imax,7), wrk2d(len), bcs(len,2)
  TREAL error, dummy
  
  TINTEGER type

! ###################################################################
  scalex = C_1_R
!  scalex = M_REAL(imax-1)
  jmax = 1
  kmax = 1
  imax_total = imax
  imode_fdm = 6
  iunif = 1
  inb_scal = 0

!  WRITE(*,*) 'Periodic (0) or nonperiodic (1) case ?'
!  READ(*,*) i1bc
  i1bc  = 1

  WRITE(*,*) 'Parameter ?'
  READ(*,*) lambda
  
  IF ( i1bc .EQ. 0 ) THEN
     DO i = 1,imax
        x(i) = M_REAL(i-1)/M_REAL(imax)*scalex
     ENDDO

  ELSE
     DO i = 1,imax
        x(i) = M_REAL(i-1)/M_REAL(imax-1)*scalex
     ENDDO
  ENDIF
  ! OPEN(21,file='grid.dat')
  ! DO i = 1,imax
  !    READ(21,*) x(i)
  ! ENDDO
  ! scalex=x(imax)-x(1)
  ! CLOSE(21)
  
  CALL DNS_INITIALIZE
  ! TO BE REVIEWED
  ! CALL FDM_INITIALIZE(iunif, imode_fdm, imax, i1bc, scalex, x, dx, wrk1d)
  x_0 = C_05_R
      
! ###################################################################
! Define the function
  DO i = 1,imax
     DO l = 1,len
! single-mode
        ! u(l,i)     =                                 &
        !       SIN(C_2_R*C_PI_R/scalex*lambda*x(i)+C_PI_R/C_4_R)
        ! du1_a(l,i) = (C_2_R*C_PI_R/scalex*lambda)    &
        !      *COS(C_2_R*C_PI_R/scalex*lambda*x(i)+C_PI_R/C_4_R)
        ! du2_a(l,i) =-(C_2_R*C_PI_R/scalex*lambda)**2 &
        !      *u(l,i)
! Gaussian
        dummy = C_1_R / ( C_2_R*(scalex/M_REAL(lambda*l))**2 )
        u(l,i)     = EXP(-dummy*(x(i)-x_0*scalex)**2)
        du1_a(l,i) =-C_2_R *dummy *(x(i)-x_0*scalex) *u(l,i)
        du2_a(l,i) =-C_2_R *dummy *(x(i)-x_0*scalex) *du1_a(l,i) - C_2_R *dummy *u(l,i)
! Exponential
        ! u(l,i)     = EXP(x(i)/(scalex/lambda))
        ! du1_a(l,i) = lambda/scalex*u(l,i)
        ! du2_a(l,i) = lambda/scalex*du1_a(l,i)
! delta-function
!  u(i)     = MAX(C_0_R,C_2_R-M_REAL(i))
!  du1_a(i) = C_0_R
!  du2_a(i) = C_0_R
! hyperboic tangent
        ! u(l,i)     = lambda*LOG(C_1_R+EXP(x(i)/lambda))
        ! du1_a(l,i) = C_05_R*(C_1_R+TANH(C_05_R*x(i)/lambda))
        ! du2_a(l,i) = C_025_R/lambda/(COSH(C_05_R*x(i)/lambda))**2
! Polynomial
        ! dummy = C_4_R
        ! u(l,i)     =                       ( (scalex-x(i)) /lambda)** dummy
        ! du1_a(l,i) = dummy                *( (scalex-x(i)) /lambda)**(dummy-C_1_R)
        ! du2_a(l,i) = dummy *(dummy-C_1_R) *( (scalex-x(i)) /lambda)**(dummy-C_2_R)
     ENDDO
  ENDDO

! ###################################################################
  type = 1

! ###################################################################
  IF ( type .EQ. 1 ) THEN ! Testing second-order derivatives
! Jacobian based
     CALL PARTIAL_XX(i1, iunif,imode_fdm, imax,jmax,kmax, i1bc,&
          dx, u, du2_n2, i0,i0, i0,i0, du1_n, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, du1_n, du2_n1, i0,i0, wrk1d,wrk2d,wrk3d)
     
! Direct metrics
     CALL FDM_C2N6N_INITIALIZE(imax, x, wrk1d(1,1), wrk1d(1,4))
     CALL TRIDFS(imax,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))

     CALL FDM_C2N6N_RHS(imax,len, wrk1d(1,4), u, du2_n3)
     CALL TRIDSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du2_n3)

  ENDIF

! -------------------------------------------------------------------
  IF ( type .EQ. 2 ) THEN ! Testing new BCs routines
     ibc = 0
     CALL FDM_C1N6_BCS_LHS(imax,     ibc, dx, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
     CALL FDM_C1N6_BCS_RHS(imax,len, ibc,     u,du1_n)
     
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
        ! WRITE(20,1000) x(i), u(l,i), du1_a(l,i), du1_n(l,i), du1_a(l,i)-du1_n(l,i)
        ! du1_n(l,i)= ABS(du1_a(l,i)-du1_n(l,i))
        ! dummy   = dummy   + du1_a(l,i)*du1_a(l,i)
        WRITE(20,1000) x(i), u(l,i), du2_a(l,i), du2_n1(l,i), du2_a(l,i)-du2_n1(l,i)
        du1_n(l,i)= ABS(du2_a(l,i)-du2_n1(l,i))
        dummy = dummy + du2_a(l,i)*du2_a(l,i)

        error = error + du1_n(l,i)*du1_n(l,i)
     ENDDO
  ENDDO
  CLOSE(20)
  
  WRITE(*,*) 'Solution L2-norm ...........:', SQRT(dx(1,1)*dummy) / M_REAL(len) 
  IF ( dummy .EQ. C_0_R ) STOP

  WRITE(*,*) 'Relative Error L2-norm .....:', SQRT(dx(1,1)*error)     /MAXVAL(ABS(du2_a))
  WRITE(*,*) 'Relative Error Linf-norm ...:', MAXVAL(du1_n(1,1:imax)) /MAXVAL(ABS(du2_a))

  STOP

1000 FORMAT(5(1x,e12.5)) 

END PROGRAM VPARTIAL
