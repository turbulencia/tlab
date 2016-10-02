#include "types.h"
#include "dns_const.h"

!# Compute dx/di and create LU factorization for first- and second-order derivatives

SUBROUTINE FDM_INITIALIZE(imethod, x, g, wrk1d)

  USE DNS_TYPES,  ONLY : grid_structure
  USE DNS_GLOBAL, ONLY : inb_scal
  USE DNS_GLOBAL, ONLY : reynolds, schmidt

  IMPLICIT NONE
  
#include "integers.h"
  
  TYPE(grid_structure),                        INTENT(INOUT) :: g
  TINTEGER,                                    INTENT(IN)    :: imethod
  TREAL, DIMENSION(g%size,g%inb_grid), TARGET, INTENT(INOUT) :: x
  TREAL, DIMENSION(g%size,5),                  INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------
  TINTEGER i, ip, is, ig, ibc_min, ibc_max, nx
  TREAL r04, r28, r24, r48, r25, r60, dummy

! ###################################################################
  nx = g%size
  
  ig = 2 ! Accumulating counter to define pointers inside array x
  
! ###################################################################
! Jacobians
! ###################################################################
  g%jac => x(:,ig:)

  IF ( nx .EQ. 1 ) THEN
     g%jac(1,1) = C_1_R
     g%jac(1,2) = C_1_R
     RETURN
  ENDIF

! ###################################################################
  IF ( g%uniform ) THEN
! -------------------------------------------------------------------
! first derivative
! -------------------------------------------------------------------
     DO i = 2,nx-1
        g%jac(i,1) = (x(i+1,1)-x(i-1,1))*C_05_R
     ENDDO

! Boundary points
     IF ( g%periodic ) THEN
        g%jac(nx,1) = ( x(1,1) + g%scale - x(nx-1,1)           )*C_05_R
        g%jac(1, 1) = ( x(2,1)-x(1,1) + x(1,1)+g%scale-x(nx,1) )*C_05_R

     ELSE
        g%jac(1, 1) = g%jac(2,1)
        g%jac(nx,1) = g%jac(nx-1,1)

     ENDIF

! -------------------------------------------------------------------
! second derivative is zero
! -------------------------------------------------------------------
     g%jac(:,2) = C_0_R

! ###################################################################
  ELSE ! derivative wrt computational grid, uniform
! -------------------------------------------------------------------
! first derivative
! -------------------------------------------------------------------
     g%jac(:,1) = C_1_R

     SELECT CASE( imethod )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4_LHS(nx,    i0,i0, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N4_RHS(nx,i1, i0,i0, x, g%jac(1,1))
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT )
        CALL FDM_C1N6_LHS(nx,    i0,i0, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N6_RHS(nx,i1, i0,i0, x, g%jac(1,1))

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8_LHS(nx,    i0,i0, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N8_RHS(nx,i1, i0,i0, x, g%jac(1,1))
        
     END SELECT

     CALL TRIDFS(nx,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
     CALL TRIDSS(nx, i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), g%jac(1,1))

! -------------------------------------------------------------------
! second derivative
! -------------------------------------------------------------------
     wrk1d(:,4) = C_1_R; wrk1d(:,5) = C_0_R

     SELECT CASE( imethod )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4_LHS(        nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N4_RHS(.TRUE., nx,i1, i0,i0, wrk1d(1,4), x,wrk1d(1,5),g%jac(1,2))
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT )
        CALL FDM_C2N6_LHS(        nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N6_RHS(.TRUE., nx,i1, i0,i0, wrk1d(1,4), x,wrk1d(1,5),g%jac(1,2))
        
     CASE( FDM_COM8_JACOBIAN ) ! Not yet developed; default to 6. order
        CALL FDM_C2N6_LHS(        nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N6_RHS(.TRUE., nx,i1, i0,i0, wrk1d(1,4), x,wrk1d(1,5),g%jac(1,2))
        
     END SELECT

     CALL TRIDFS(nx,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))
     CALL TRIDSS(nx, i1, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), g%jac(1,2))

  ENDIF

  ig = ig + 2

! ###################################################################
! LU factorization first-order derivative, done in routine TRID*FS
! ###################################################################
  g%lu1 => x(:,ig:)

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
  IF ( g%periodic ) THEN
     SELECT CASE( imethod )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4P_LHS(nx, g%jac, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3))
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C1N6P_LHS(nx, g%jac, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3))
        
     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8P_LHS(nx, g%jac, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3))
        
     END SELECT
     
     CALL TRIDPFS(nx, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5))
     ig = ig + 5

! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal for 4 different BCs
! -------------------------------------------------------------------
  ELSE
     DO i = 0,3
        ibc_min = MOD(i,2)
        ibc_max = i /2
        ip = i*3

        SELECT CASE( imethod )
           
        CASE( FDM_COM4_JACOBIAN )
           CALL FDM_C1N4_LHS(nx, ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        CASE( FDM_COM6_JACOBIAN )
           CALL FDM_C1N6_LHS(nx, ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        CASE( FDM_COM8_JACOBIAN )
           CALL FDM_C1N8_LHS(nx, ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        CASE( FDM_COM6_DIRECT   ) ! Not yet implemented; using Jacobian version
           CALL FDM_C1N6_LHS(nx, ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        END SELECT
        
        CALL TRIDFS(nx, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))
        ig = ig + 3
     ENDDO
  
  ENDIF

! ###################################################################
! LU factorization second-order derivative, done in routine TRID*FS
! ###################################################################
  g%lu2 => x(:,ig:)

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
  IF ( g%periodic ) THEN
     SELECT CASE( imethod )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4P_LHS(nx, g%jac, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3))
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT )  ! Direct = Jacobian because uniform grid
        CALL FDM_C2N6P_LHS(nx, g%jac, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3))

     CASE( FDM_COM8_JACOBIAN )                   ! Not yet developed
        CALL FDM_C2N6P_LHS(nx, g%jac, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3))

     END SELECT
     
     CALL TRIDPFS(nx, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3),g%lu2(1,4),g%lu2(1,5))
     ig = ig + 5
     
! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal for 4 different BCs
! -------------------------------------------------------------------
  ELSE
     DO i = 0,3
        ibc_min = MOD(i,2)
        ibc_max = i /2
        ip = i*3
        SELECT CASE( imethod )
           
        CASE( FDM_COM4_JACOBIAN )
           CALL FDM_C2N4_LHS(nx, ibc_min,ibc_max, g%jac, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3))

        CASE( FDM_COM6_JACOBIAN )
           
           CALL FDM_C2N6_LHS(nx, ibc_min,ibc_max, g%jac, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3))

        CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented
           CALL FDM_C2N6_LHS(nx, ibc_min,ibc_max, g%jac, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3)) ! 8th not yet developed

        CASE( FDM_COM6_DIRECT   )
           IF ( i .EQ. 0 ) CALL FDM_C2N6N_INITIALIZE(nx, x, g%lu2(1,ip+1), g%lu2(1,ip+4))
           
        END SELECT

! The direct mode is only implemented for bcs=(0,0); we use the remaining array
! to save other data        
        IF  ( imethod .EQ. FDM_COM6_DIRECT ) THEN
           IF ( i .EQ. 0 ) THEN
              g%lu2(:,ip+8:ip+10) = g%lu2(:,ip+1:ip+3) ! saving the array A w/o LU decomposition
              CALL TRIDFS(nx, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3))
              ig = ig + 10
           ENDIF
        ELSE
           CALL TRIDFS(nx, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3))
           ig = ig + 3
        ENDIF

     ENDDO

  ENDIF

! ###################################################################
! LU factorization second-order derivative times the diffusivities
! ###################################################################
  g%lu2d => x(:,ig:)

  ip = 0
  DO is = 0,inb_scal ! case 0 for the reynolds number
     IF ( is .EQ. 0 ) THEN; dummy = reynolds
     ELSE;                  dummy = reynolds*schmidt(is); ENDIF

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
     IF ( g%periodic ) THEN ! Check routines TRIDPFS and TRIDPSS
        g%lu2d(:,ip+1) = g%lu2(:,1)          ! matrix L; 1. subdiagonal
        g%lu2d(:,ip+2) = g%lu2(:,2) /dummy   ! matrix L; 1/diagonal
        g%lu2d(:,ip+3) = g%lu2(:,3)          ! matrix U is the same
        g%lu2d(:,ip+4) = g%lu2(:,4) *dummy   ! matrix L; Additional row/column         
        g%lu2d(:,ip+5) = g%lu2(:,5)          ! matrix U is the same
           
        ig = ig + 5
        ip = ip + 5

! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal, 1 single BCs
! -------------------------------------------------------------------
     ELSE                   ! Check routines TRIDFS and TRIDSS
        g%lu2d(:,ip+1) = g%lu2(:,1)          ! matrix L is the same
        g%lu2d(:,ip+2) = g%lu2(:,2) /dummy   ! matrix U; 1/diagonal
        g%lu2d(:,ip+3) = g%lu2(:,3) *dummy   ! matrix U; 1. superdiagonal

        ig = ig + 3
        ip = ip + 3
        
     ENDIF
     
  ENDDO

! ###################################################################
! Modified wavenumbers in periodic case
! ###################################################################
  g%mwn => x(:,ig:)

  IF ( g%periodic ) THEN
     DO i = 1,nx ! Define wavenumbers
        IF ( i .LE. nx/2+1 ) THEN; wrk1d(i,1) = C_2_R*C_PI_R*M_REAL(i-1   )/M_REAL(nx)
        ELSE;                      wrk1d(i,1) = C_2_R*C_PI_R*M_REAL(i-1-nx)/M_REAL(nx); ENDIF
     ENDDO
        
! -------------------------------------------------------------------
! First-order derivative
! -------------------------------------------------------------------
     r04 = C_2_R /C_5_R
     r28 = C_14_R*C_2_R
     r48 = C_6_R *C_8_R
     r25 = C_5_R *C_5_R /C_4_R
     r60 = C_1_R /( C_6_R *C_10_R )

     SELECT CASE( imethod )

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT )
        g%mwn(:,1)=( r28*sin(wrk1d(:,1))+    sin(C_2_R*wrk1d(:,1))                          )&
                  /( C_18_R +C_12_R*cos(wrk1d(:,1)))

     CASE( FDM_COM8_JACOBIAN )
        g%mwn(:,1)=( r25*sin(wrk1d(:,1))+r04*sin(C_2_R*wrk1d(:,1))-r60*sin(C_3_R*wrk1d(:,1)))&
                  /( C_4_R +C_3_R *cos(wrk1d(:,1)))

     END SELECT

! Final calculations because it is mainly used in the Poisson solver like this
     g%mwn(:,1) = ( g%mwn(:,1) /g%jac(1,1) )**2

! -------------------------------------------------------------------
! Second-order derivative
! -------------------------------------------------------------------
     r24 = C_6_R*C_4_R 

     SELECT CASE( imethod )

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT )
        g%mwn(:,2)=( r24*(1-cos(wrk1d(:,1))) + C_1_5_R*(C_1_R-cos(C_2_R*wrk1d(:,1))) )&
                  /( C_11_R + C_4_R*cos(wrk1d(:,1)) )

     CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented

     END SELECT

! Final calculations because it is mainly used in the Helmholtz solver like this
     g%mwn(:,2) = g%mwn(:,2) /( g%jac(1,1)**2 )

     ig = ig + 2

  ENDIF

! ###################################################################
! Check array sizes
! ###################################################################
  ! IF ( ig .NE. g%inb_grid ) THEN
  !    CALL IO_WRITE_ASCII(efile, 'FDM_INITIALIZE. Grid size incorrect.')
  !    CALL DNS_STOP(DNS_ERROR_DIMGRID)
  ! ENDIF
  
  RETURN
END SUBROUTINE FDM_INITIALIZE
