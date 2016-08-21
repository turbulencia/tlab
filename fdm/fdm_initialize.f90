#include "types.h"
#include "dns_const.h"

!# Compute dx/di and create LU factorization for first- and second-order derivatives

SUBROUTINE FDM_INITIALIZE(imethod, nx, scalex, uniform, periodic, x, dx, wrk1d)

  USE DNS_GLOBAL, ONLY : inb_scal
  USE DNS_GLOBAL, ONLY : inb_grid, inb_grid_1, inb_grid_2, inb_grid_3
  USE DNS_GLOBAL, ONLY : reynolds, schmidt

  IMPLICIT NONE
  
#include "integers.h"
  
  TINTEGER,                    INTENT(IN) :: imethod, nx
  TREAL,                       INTENT(IN) :: scalex
  LOGICAL,                     INTENT(IN) :: uniform, periodic
  TREAL,DIMENSION(nx),         INTENT(IN) :: x
  TREAL,DIMENSION(nx,inb_grid),INTENT(OUT):: dx
  TREAL,DIMENSION(nx,5)                   :: wrk1d

! -------------------------------------------------------------------
  TINTEGER i, ip, is
  TREAL kx, r28, r24, r48, r25, r60, dummy

! ###################################################################
  IF ( nx .EQ. 1 ) THEN
     dx(1,1) = C_1_R
     dx(1,2) = C_1_R
     RETURN
  ENDIF

! ###################################################################
! Uniform grid
! ###################################################################
  IF ( uniform ) THEN
! -------------------------------------------------------------------
! first derivative
! -------------------------------------------------------------------
     DO i = 2,nx-1
        dx(i,1) = (x(i+1)-x(i-1))*C_05_R
     ENDDO
! periodic BCs
     IF ( periodic ) THEN
        dx(nx,1) = ( x(1) + scalex - x(nx-1)       )*C_05_R
        dx(1, 1) = ( x(2)-x(1) + x(1)+scalex-x(nx) )*C_05_R
! nonperiodic BCs
     ELSE
        dx(1, 1) = dx(2,1)
        dx(nx,1) = dx(nx-1,1)
     ENDIF

! -------------------------------------------------------------------
! second derivative is zero
! -------------------------------------------------------------------
     dx(:,2) = C_0_R

! ###################################################################
! Nonuniform grid
! ###################################################################
  ELSE
! -------------------------------------------------------------------
! first derivative
! -------------------------------------------------------------------
     dx(:,1) = C_1_R

     IF      ( imethod .EQ. FDM_COM4_JACOBIAN                                 ) THEN
        CALL FDM_C1N4_LHS(nx,    i0,i0, dx, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N4_RHS(nx,i1, i0,i0, x, dx(1,1))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN .OR. imethod .EQ. FDM_COM6_DIRECT ) THEN
        CALL FDM_C1N6_LHS(nx,    i0,i0, dx, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N6_RHS(nx,i1, i0,i0, x, dx(1,1))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN                                 ) THEN
        CALL FDM_C1N8_LHS(nx,    i0,i0, dx, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C1N8_RHS(nx,i1, i0,i0, x, dx(1,1))
     ENDIF

     CALL TRIDFS(nx,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
     CALL TRIDSS(nx, i1, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), dx(1,1))

! -------------------------------------------------------------------
! second derivative
! -------------------------------------------------------------------
     wrk1d(:,4) = C_1_R; wrk1d(:,5) = C_0_R

! derivative wrt computational grid, uniform
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN                                 ) THEN
        CALL FDM_C2N4_LHS(    nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N4_RHS(i0, nx,i1, i0,i0, wrk1d(1,4), x,wrk1d(1,5),dx(1,2))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN .OR. imethod .EQ. FDM_COM6_DIRECT ) THEN
        CALL FDM_C2N6_LHS(    nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N6_RHS(i0, nx,i1, i0,i0, wrk1d(1,4), x,wrk1d(1,5),dx(1,2))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN                                 ) THEN !8th not yet developed
        CALL FDM_C2N6_LHS(    nx,    i0,i0, wrk1d(1,4), wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL FDM_C2N6_RHS(i0, nx,i1, i0,i0, wrk1d(1,4), x,wrk1d(1,5),dx(1,2))
     ENDIF

     CALL TRIDFS(nx,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))
     CALL TRIDSS(nx, i1, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), dx(1,2))

  ENDIF

! ###################################################################
! LU factorization first-order derivative, done in routine TRID*FS
! ###################################################################
  ip = inb_grid_1 - 1

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
  IF ( periodic ) THEN
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN                                 ) THEN
        CALL FDM_C1N4P_LHS(nx, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN .OR. imethod .EQ. FDM_COM6_DIRECT ) THEN
        CALL FDM_C1N6P_LHS(nx, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN                                 ) THEN
        CALL FDM_C1N8P_LHS(nx, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ENDIF

     CALL TRIDPFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5))

! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal for 4 different BCs
! -------------------------------------------------------------------
  ELSE
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C1N4_LHS(nx, i0,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C1N6_LHS(nx, i0,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C1N8_LHS(nx, i0,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN ! not yet implemented
        CALL FDM_C1N6_LHS(nx, i0,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

     ip = ip + 3
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C1N4_LHS(nx, i1,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C1N6_LHS(nx, i1,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C1N8_LHS(nx, i1,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN ! not yet implemented
        CALL FDM_C1N6_LHS(nx, i1,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

     ip = ip + 3
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C1N4_LHS(nx, i0,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C1N6_LHS(nx, i0,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C1N8_LHS(nx, i0,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN ! not yet implemented
        CALL FDM_C1N6_LHS(nx, i0,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

     ip = ip + 3
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C1N4_LHS(nx, i1,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C1N6_LHS(nx, i1,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C1N8_LHS(nx, i1,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN ! not yet implemented
        CALL FDM_C1N6_LHS(nx, i1,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

  ENDIF

! ###################################################################
! LU factorization second-order derivative, done in routine TRID*FS
! ###################################################################
  ip = inb_grid_2 - 1

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
  IF ( periodic ) THEN
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN                                 ) THEN
        CALL FDM_C2N4P_LHS(nx, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN .OR. imethod .EQ. FDM_COM6_DIRECT ) THEN
        CALL FDM_C2N6P_LHS(nx, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN                                 ) THEN
        CALL FDM_C2N6P_LHS(nx, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3)) ! 8th not yet developed
     ENDIF
     CALL TRIDPFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5))

! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal
! -------------------------------------------------------------------
  ELSE
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C2N4_LHS(nx, i0,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i0,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i0,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3)) ! 8th not yet developed
     ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN
        CALL FDM_C2N6N_INITIALIZE(nx, x, dx(1,ip+1), dx(1,ip+4))
        dx(:,ip+8:ip+10) = dx(:,ip+1:ip+3) ! saving the array A w/o LU decomposition
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

! Different BCs
     IF ( imethod .NE. FDM_COM6_DIRECT ) THEN

     ip = ip + 3
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C2N4_LHS(nx, i1,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i1,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i1,i0, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3)) ! 8th not yet developed
     ! ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN ! default to standard case
     !    dx(1,ip+1:ip+3) = dx(1,inb_grid_2+7:inb_grid_2+9)
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

     ip = ip + 3
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C2N4_LHS(nx, i0,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i0,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i0,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3)) ! 8th not yet developed
     ! ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN ! default to standard case
     !    dx(1,ip+1:ip+3) = dx(1,inb_grid_2+7:inb_grid_2+9)
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

     ip = ip + 3
     IF      ( imethod .EQ. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C2N4_LHS(nx, i1,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i1,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))
     ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C2N6_LHS(nx, i1,i1, dx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3)) ! 8th not yet developed
     ! ELSE IF ( imethod .EQ. FDM_COM6_DIRECT   ) THEN ! default to standard case
     !    dx(1,ip+1:ip+3) = dx(1,inb_grid_2+7:inb_grid_2+9)
     ENDIF
     CALL TRIDFS(nx, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3))

  ENDIF

  ENDIF

! ###################################################################
! LU factorization second-order derivative times the diffusivities
! ###################################################################
  ip = inb_grid_3 - 1

  DO is = 0,inb_scal ! case 0 for the reynolds number
     IF ( is .EQ. 0 ) THEN; dummy = reynolds
     ELSE;                  dummy = reynolds*schmidt(is); ENDIF

! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
     IF ( periodic ) THEN ! Check routines TRIDPFS and TRIDPSS
        dx(:,ip+1) = dx(:,inb_grid_2  )         ! matrix L; 1. subdiagonal
        dx(:,ip+2) = dx(:,inb_grid_2+1) /dummy  ! matrix L; 1/diagonal
        dx(:,ip+3) = dx(:,inb_grid_2+2)         ! matrix U is the same
        dx(:,ip+4) = dx(:,inb_grid_2+3) *dummy  ! matrix L; Additional row/column         
        dx(:,ip+5) = dx(:,inb_grid_2+4)         ! matrix U is the same
           
! -------------------------------------------------------------------
! Nonperiodic case; tridiagonal, 1 single BCs
! -------------------------------------------------------------------
     ELSE                    ! Check routines TRIDFS and TRIDSS
        dx(:,ip+1) = dx(:,inb_grid_2  )         ! matrix L is the same
        dx(:,ip+2) = dx(:,inb_grid_2+1) /dummy  ! matrix U; 1/diagonal
        dx(:,ip+3) = dx(:,inb_grid_2+2) *dummy  ! matrix U; 1. superdiagonal
     ENDIF
     
     ip = ip + 5 ! stride is 5, which is set by the maximum imposed by periodic BCs
  ENDDO

! ###################################################################
! Modified wavenumbers in periodic case
! ###################################################################
  IF ( periodic ) THEN
     
! -------------------------------------------------------------------
! First order derivative
! -------------------------------------------------------------------
     ip = inb_grid_1 + 5

     r28 = C_14_R*C_2_R
     r48 = C_6_R*C_8_R
     r25 = C_5_R*C_5_R/C_4_R
     r60 = C_1_R/(C_6_R*C_10_R)
     DO i = 1,nx
        IF ( i .LE. nx/2+1 ) THEN; kx = C_2_R*C_PI_R*M_REAL(i-1   )/M_REAL(nx)
        ELSE;                      kx = C_2_R*C_PI_R*M_REAL(i-1-nx)/M_REAL(nx); ENDIF
        IF      ( imethod .EQ. FDM_COM6_JACOBIAN .OR. imethod .EQ. FDM_COM6_DIRECT ) THEN
           dx(i,ip)=(r28*sin(kx)+            sin(C_2_R*kx)                  )/(C_18_R+C_12_R*cos(kx))
        ELSE IF ( imethod .EQ. FDM_COM8_JACOBIAN                                 ) THEN
           dx(i,ip)=(r25*sin(kx)+C_2_R/C_5_R*sin(C_2_R*kx)-r60*sin(C_3_R*kx))/(C_4_R +C_3_R *cos(kx))
        ENDIF
     ENDDO

! Final calculations because it is mainly used in the Poisson solver like this
     dx(:,ip) = (dx(:,ip)/dx(1,1))**2

! -------------------------------------------------------------------
! Second order derivative
! -------------------------------------------------------------------
     ip = inb_grid_2 + 5

     r24 = C_6_R*C_4_R 
     DO i = 1,nx
        IF ( i .LE. nx/2+1 ) THEN; kx = C_2_R*C_PI_R*M_REAL(i-1   )/M_REAL(nx)
        ELSE;                      kx = C_2_R*C_PI_R*M_REAL(i-1-nx)/M_REAL(nx); ENDIF
        IF      ( imethod .EQ. FDM_COM6_JACOBIAN .OR. imethod .EQ. FDM_COM6_DIRECT ) THEN
           dx(i,ip)=( r24*(1-cos(kx)) + C_1_5_R*(C_1_R-cos(C_2_R*kx)) )/( C_11_R + C_4_R*cos(kx) )
        ENDIF
     ENDDO

! Final calculations because it is mainly used in the Helmholtz solver like this
     dx(:,ip) = dx(:,ip)/(dx(1,1)**2)

  ENDIF

  RETURN
END SUBROUTINE FDM_INITIALIZE
