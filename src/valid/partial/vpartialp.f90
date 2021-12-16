!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/16 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate  6th-order pentadiagonal compact scheme.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

PROGRAM VPARTIALP

  USE TLAB_TYPES, ONLY : grid_dt
  
  IMPLICIT NONE

#include "integers.h"
 
TYPE(grid_dt)                        :: g
  
  TINTEGER                           :: jmax,kmax, i, l
  TINTEGER, PARAMETER                :: imax=128, len=1, inb_grid=42
  
  TREAL,    DIMENSION(imax,inb_grid) :: x
  TREAL,    DIMENSION(len,imax)      :: u, wrk3d
  TREAL,    DIMENSION(len,imax)      :: du1_a, du1_b
  TREAL,    DIMENSION(len,imax)      :: du2_a
  TREAL,    DIMENSION(imax,7)        :: wrk1d
  TREAL,    DIMENSION(len)           :: wrk2d
  TREAL,    DIMENSION(len,2)         :: bcs
  
  TREAL                              :: lambda, error, dummy
  TINTEGER                           :: test_type, ibc
 
! ###################################################################
! Initialize grid
  g%size     = imax 
  g%scale    = C_1_R
  g%uniform  = .TRUE.
  jmax       = 1
  kmax       = 1
  g%periodic = .FALSE.

! Valid stettings
  lambda     = 1 ! periodicity
  test_type  = 2 
  ibc        = 3

! Scheme
  g%mode_fdm = FDM_COM6_JACPENTA 

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

! ###################################################################
! Define the function and analytic derivatives
  DO i = 1,imax
    DO l = 1,len
      u(l,i)     =                                  &
             SIN(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))
      du1_a(l,i) = (C_2_R*C_PI_R/g%scale*lambda)    &
            *COS(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))
      du2_a(l,i) =-(C_2_R*C_PI_R/g%scale*lambda)**2 &
            *u(l,i)
    ENDDO
  ENDDO

! ###################################################################
! Testing first-order derivatives
  
  IF ( test_type .EQ. 1 ) THEN 
    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g, u, du1_b, wrk3d, wrk2d,wrk3d)

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
      WRITE(20,1000) g%nodes(i), u(l,i), du1_a(l,i), du1_b(l,i), du1_a(l,i) - du1_b(l,i)
      du1_c(l,i) = ABS(du1_a(l,i) - du1_b(l,i))
      dummy = dummy + du1_a(l,i) * du1_a(l,i)
      error = error + du1_c(l,i) * du1_c(l,i)
    ENDDO
  ENDDO
  CLOSE(20)
  
  WRITE(*,*) 'Solution L2-norm ...........:', SQRT(g%jac(1,1)*dummy) / M_REAL(len)
  IF ( dummy .EQ. C_0_R ) STOP
  WRITE(*,*) 'Relative Error L2-norm .....:', SQRT(g%jac(1,1)*error)  /MAXVAL(ABS(du1_a))
  WRITE(*,*) 'Relative Error Linf-norm ...:', MAXVAL(du1_c(1,1:imax)) /MAXVAL(ABS(du1_a))
  
  STOP
  
  1000 FORMAT(5(1x,e12.5))
  
END PROGRAM VPARTIALP