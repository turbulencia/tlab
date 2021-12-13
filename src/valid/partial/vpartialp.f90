!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/13 - J. Kostelecky
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
  TREAL,    DIMENSION(len,imax)      :: du1_a, du1_b, du1_c, du1_d, du1_e
  TREAL,    DIMENSION(len,imax)      :: du2_a
  TREAL,    DIMENSION(imax,7)        :: wrk1d
  TREAL,    DIMENSION(len)           :: wrk2d
  TREAL,    DIMENSION(len,2)         :: bcs
  
  TREAL                              :: lambda, error, dummy
 
! ###################################################################
! Initialize grid
  g%size = imax; g%scale = C_1_R; g%uniform  = .TRUE.
  jmax   = 1; kmax = 1
  lambda = 1 ! periodicity

  g%mode_fdm = FDM_COM6_JACPENTA ! FDM_COM6_JACOBIAN
  g%periodic = .FALSE. ! .TRUE.
  
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
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g, u, du1_b, wrk3d, wrk2d,wrk3d)
  
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