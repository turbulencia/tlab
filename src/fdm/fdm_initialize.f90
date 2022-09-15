#include "types.h"
#include "dns_const.h"

!# Compute dx/di and create LU factorization for first- and second-order derivatives

SUBROUTINE FDM_INITIALIZE(x, g, wrk1d)
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
#endif
  USE TLAB_TYPES,  ONLY : grid_dt
  USE TLAB_VARS, ONLY : inb_scal, istagger
  USE TLAB_VARS, ONLY : reynolds, schmidt
  USE TLAB_VARS, ONLY : C1N6M_ALPHA2, C1N6M_BETA2
  USE TLAB_VARS, ONLY : C1N6M_A, C1N6M_BD2, C1N6M_CD3

  USE TLAB_PROCS

  IMPLICIT NONE

#include "integers.h"

  TYPE(grid_dt),                               INTENT(INOUT) :: g
  TREAL, DIMENSION(g%size,g%inb_grid), TARGET, INTENT(INOUT) :: x
  TREAL, DIMENSION(g%size,5),                  INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------
  TINTEGER i, ip, is, ig, ibc_min, ibc_max, nx
  TREAL r04, r28, r24, r48, r25, r60, dummy
  TREAL ra1, rb2, rc3, r2a, r2b 
  TREAL sa,  sb,  sal

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'Entering SUBROUTINE FDM_INITIALIZE')
#endif

! ###################################################################
  nx = g%size

  ig = 1 ! Accumulating counter to define pointers inside array x

! ###################################################################
  g%nodes => x(:,ig)

  ig = ig + 1

! ###################################################################
! Coefficients of pentadiagonal 6th-order scheme for 1st derivative
  IF (g%mode_fdm .EQ. FDM_COM6_JACPENTA) CALL FDM_C1N6M_COEFF()
  
! ###################################################################
! Jacobians
! ###################################################################
  g%jac => x(:,ig:)

  IF ( nx .EQ. 1 ) THEN
     g%jac(1,1) = C_1_R
     g%jac(1,2) = C_1_R
     g%jac(1,3) = C_0_R
     g%jac(1,4) = C_0_R
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

     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4_LHS(nx,    i0,i0, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N4_RHS(nx,i1, i0,i0, x, g%jac(1,1))

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT )
        CALL FDM_C1N6_LHS(nx,    i0,i0, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N6_RHS(nx,i1, i0,i0, x, g%jac(1,1))
      
     CASE( FDM_COM6_JACPENTA )
        CALL FDM_C1N6M_LHS(nx,    i0,i0, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL FDM_C1N6M_RHS(nx,i1, i0,i0, x, g%jac(1,1))

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8_LHS(nx,    i0,i0, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N8_RHS(nx,i1, i0,i0, x, g%jac(1,1))

     END SELECT

     IF (.NOT. (g%mode_fdm .EQ. FDM_COM6_JACPENTA)) THEN
        CALL TRIDFS(nx,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(nx, i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), g%jac(1,1))
     ELSE
        CALL PENTADFS2(nx,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL PENTADSS2(nx, i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), g%jac(1,1))
     ENDIF

! -------------------------------------------------------------------
! second derivative
! -------------------------------------------------------------------
     wrk1d(:,4) = C_1_R; wrk1d(:,5) = C_0_R

     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4_LHS(nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N4_RHS(nx,i1, i0,i0,             x, g%jac(1,2))

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA )
       ! CALL FDM_C2N6_LHS( nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
       ! CALL FDM_C2N6_RHS( nx,i1, i0,i0,             x, g%jac(1,2))
       CALL FDM_C2N6H_LHS( nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
       CALL FDM_C2N6H_RHS( nx,i1, i0,i0,             x, g%jac(1,2))

     CASE( FDM_COM8_JACOBIAN ) ! Not yet developed; default to 6. order
        CALL FDM_C2N6_LHS( nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N6_RHS( nx,i1, i0,i0,             x, g%jac(1,2))

     END SELECT

     CALL TRIDFS(nx,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))
     CALL TRIDSS(nx, i1, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), g%jac(1,2))

  ENDIF

! ###################################################################
! Saving operations for the time-stability constraint
  g%jac(:,3) = C_1_R /g%jac(:,1)
  g%jac(:,4) = g%jac(:,3) *g%jac(:,3)

  ig = ig + 4

! ###################################################################
! LU factorization first-order derivative, done in routine TRID*FS
! ###################################################################
  g%lu1 => x(:,ig:)

! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( g%periodic ) THEN
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4P_LHS(nx, g%jac, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3))

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C1N6P_LHS(nx, g%jac, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3))

     CASE( FDM_COM6_JACPENTA )
        CALL FDM_C1N6MP_LHS(nx, g%jac, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5))     

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8P_LHS(nx, g%jac, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3))

     END SELECT

     IF (.NOT. (g%mode_fdm .EQ. FDM_COM6_JACPENTA)) THEN
        CALL TRIDPFS(  nx, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5))
     ELSE
        CALL PENTADPFS(nx, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5),g%lu1(1,6),g%lu1(1,7))
     ENDIF
     ig = ig + 7

! -------------------------------------------------------------------
! Nonperiodic case (4 different BCs)
! -------------------------------------------------------------------
  ELSE
     DO i = 0,3
        ibc_min = MOD(i,2)
        ibc_max = i /2
        ip = i*5
        SELECT CASE( g%mode_fdm )

        CASE( FDM_COM4_JACOBIAN )
           CALL FDM_C1N4_LHS(nx,  ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        CASE( FDM_COM6_JACOBIAN )
           CALL FDM_C1N6_LHS(nx,  ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        CASE( FDM_COM6_JACPENTA )
           CALL FDM_C1N6M_LHS(nx, ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3),g%lu1(1,ip+4),g%lu1(1,ip+5))

        CASE( FDM_COM8_JACOBIAN )
           CALL FDM_C1N8_LHS(nx,  ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        CASE( FDM_COM6_DIRECT   ) ! Not yet implemented; using Jacobian version
           CALL FDM_C1N6_LHS(nx,  ibc_min,ibc_max, g%jac, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))

        END SELECT

        IF (.NOT. (g%mode_fdm .EQ. FDM_COM6_JACPENTA)) THEN
           CALL TRIDFS(   nx, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3))
        ELSE
           CALL PENTADFS2(nx, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3),g%lu1(1,ip+4),g%lu1(1,ip+5))
        ENDIF
        ig = ig + 5
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
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4P_LHS(nx, g%jac, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3))

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA )  ! Direct = Jacobian because uniform grid
       ! CALL FDM_C2N6P_LHS(nx, g%jac, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3))
       CALL FDM_C2N6HP_LHS(nx, g%jac, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3))

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
        SELECT CASE( g%mode_fdm )

        CASE( FDM_COM4_JACOBIAN )
           CALL FDM_C2N4_LHS(nx, ibc_min,ibc_max, g%jac, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3))

        CASE( FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA )
          ! CALL FDM_C2N6_LHS(nx, ibc_min,ibc_max, g%jac, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3))
          CALL FDM_C2N6H_LHS(nx, ibc_min,ibc_max, g%jac, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3))

        CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented
           CALL FDM_C2N6_LHS(nx, ibc_min,ibc_max, g%jac, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3)) ! 8th not yet developed

        CASE( FDM_COM6_DIRECT   )
           IF ( i .EQ. 0 ) CALL FDM_C2N6ND_INITIALIZE(nx, x, g%lu2(1,ip+1), g%lu2(1,ip+4))

        END SELECT

! The direct mode is only implemented for bcs=(0,0); we use the remaining array
! to save other data
        IF  ( g%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
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
! LU factorization interpolation, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
  IF ( (istagger .EQ. 1) .AND. g%periodic ) THEN
     g%lu0i => x(:,ig:)

     SELECT CASE( g%mode_fdm )
     CASE DEFAULT
        CALL FDM_C0INT6P_LHS(nx, g%lu0i(1,1),g%lu0i(1,2),g%lu0i(1,3))
     END SELECT
     CALL TRIDPFS(nx, g%lu0i(1,1),g%lu0i(1,2),g%lu0i(1,3),g%lu0i(1,4),g%lu0i(1,5))
     ig = ig + 5
  ENDIF

! ###################################################################
! LU factorization first interp. derivative, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
  IF ( (istagger .EQ. 1) .AND. g%periodic ) THEN
     g%lu1i => x(:,ig:)
 
     SELECT CASE( g%mode_fdm )
     CASE DEFAULT
        CALL FDM_C1INT6P_LHS(nx, g%jac, g%lu1i(1,1),g%lu1i(1,2),g%lu1i(1,3))
     END SELECT
     CALL TRIDPFS(nx, g%lu1i(1,1),g%lu1i(1,2),g%lu1i(1,3),g%lu1i(1,4),g%lu1i(1,5))
     ig = ig + 5
  ENDIF

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
     ! tridiagonal scheme
     r04 = C_2_R /C_5_R
     r28 = C_14_R*C_2_R
     r48 = C_6_R *C_8_R
     r25 = C_5_R *C_5_R /C_4_R
     r60 = C_1_R /( C_6_R *C_10_R )
     ! pentadiagonal scheme
     ra1 = C1N6M_A
     rb2 = C1N6M_BD2
     rc3 = C1N6M_CD3
     r2a = C1N6M_ALPHA2
     r2b = C1N6M_BETA2
     ! staggered tridiagonal 6th-order scheme
     sal = C_9_R  / C_62_R 
     sa  = C_63_R / C_62_R
     sb  = C_17_R / C_62_R 

     IF ( istagger  .EQ. 0 ) THEN

        SELECT CASE( g%mode_fdm )
         
        CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT )
           g%mwn(:,1)=( r28*sin(wrk1d(:,1))+    sin(C_2_R*wrk1d(:,1))                          )&
                     /( C_18_R +C_12_R*cos(wrk1d(:,1)))
   
        CASE( FDM_COM6_JACPENTA )
           g%mwn(:,1)=( ra1*sin(wrk1d(:,1))+rb2*sin(C_2_R*wrk1d(:,1))+rc3*sin(C_3_R*wrk1d(:,1)))&
                     /( C_1_R  +r2a   *cos(wrk1d(:,1))+r2b*cos(2*wrk1d(:,1)))
         
        CASE( FDM_COM8_JACOBIAN )
           g%mwn(:,1)=( r25*sin(wrk1d(:,1))+r04*sin(C_2_R*wrk1d(:,1))-r60*sin(C_3_R*wrk1d(:,1)))&
                     /( C_4_R  +C_3_R *cos(wrk1d(:,1)))
   
        END SELECT

     ELSE ! staggered case has different modified wavenumbers!

        SELECT CASE( g%mode_fdm )
         
        CASE DEFAULT
           g%mwn(:,1)=( C_2_R * sa * sin(C_1_R/C_2_R * wrk1d(:,1)) + C_2_R/C_3_R * sb * sin(C_3_R/C_2_R * wrk1d(:,1)))&
                     /( C_1_R + C_2_R * sal * cos(wrk1d(:,1)))
 
        END SELECT

     ENDIF

! Final calculations because it is mainly used in the Poisson solver like this
     g%mwn(:,1) = ( g%mwn(:,1) /g%jac(1,1) )**2

! -------------------------------------------------------------------
! Second-order derivative
! -------------------------------------------------------------------
     r24 = C_6_R*C_4_R

     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA )
        g%mwn(:,2)=( r24*(1-cos(wrk1d(:,1))) + C_1_5_R*(C_1_R-cos(C_2_R*wrk1d(:,1))) )&
                  /( C_11_R + C_4_R*cos(wrk1d(:,1)) )

     CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented

     END SELECT

! Final calculations because it is mainly used in the Helmholtz solver like this
     g%mwn(:,2) = g%mwn(:,2) /( g%jac(1,1)**2 )

     ig = ig + 2

  ENDIF

! ###################################################################
! Density correction in anelastic mode
! ###################################################################
  g%rhoinv => x(:,ig)

  g%anelastic = .FALSE. ! Default; activated in FI_BACKGROUND_INITIALIZE

  ig = ig +1
! ###################################################################
! Check array sizes
! ###################################################################
  ! IF ( ig .NE. g%inb_grid ) THEN
  !    CALL TLAB_WRITE_ASCII(efile, 'FDM_INITIALIZE. Grid size incorrect.')
  !    CALL TLAB_STOP(DNS_ERROR_DIMGRID)
  ! ENDIF

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'Leaving SUBOURINTE FDM_INITIALIZE')
#endif

  RETURN
END SUBROUTINE FDM_INITIALIZE
