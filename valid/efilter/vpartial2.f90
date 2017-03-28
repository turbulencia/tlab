#include "types.h"

PROGRAM VPARTIAL2

  USE DNS_TYPES, ONLY : grid_dt
  IMPLICIT NONE

#include "integers.h"

  TYPE(grid_dt),            INTENT(IN)    :: g
  TINTEGER imode_fdm, imax, jmax, kmax, i, wk, i1bc, idummy, iunif
  PARAMETER(imax=128)
  TREAL scalex
  TREAL x(imax), dx(imax,2+4*3+4*3), u(imax), du1(imax), du2(imax), due(imax)
  TREAL wrk1d(imax,5), wrk2d(imax), wrk3d(imax)
  TREAL tmp(imax)

! ###################################################################
  scalex = C_2_R*C_PI_R
  jmax = 1
  kmax = 1
  imode_fdm = 6
  iunif = 1

  WRITE(*,*) 'Periodic (0) or nonperiodic (1) case ?'
  READ(*,*) i1bc
  WRITE(*,*) 'Wavenumber ?'
  READ(*,*) wk

! CHANGE TO UPDATE NEW GRID_DT
  IF ( i1bc .EQ. 0 ) THEN
     DO i = 1,imax
        x(i) = M_REAL(i-1)/M_REAL(imax)*scalex
     ENDDO
  ELSE
     OPEN(21,file='y.dat')
     DO i = 1,imax
!        x(i) = M_REAL(i-1)/M_REAL(imax-1)*scalex
        READ(21,*) idummy, x(i)
     ENDDO
     CLOSE(21)
     scalex = x(imax)-x(1)
  ENDIF

  CALL FDM_INITIALIZE(x, g, wrk1d)

! ###################################################################
! Define the function
  DO i = 1,imax
     u(i) = SIN(C_2_R*C_PI_R/scalex*M_REAL(wk)*x(i))
     due(i) = -(C_2_R*C_PI_R/scalex*M_REAL(wk))**2*u(i)
!     u(i) = EXP(-(x(i)-C_PI_R)**2/(C_PI_R**2/64.))
!     due(i) = -C_2_R*(x(i)-C_PI_R)/(C_PI_R**2/64.)*u(i)
  ENDDO

! ###################################################################
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp, du1, i1, i1, wrk1d, wrk2d, wrk3d)

  CALL PARTIAL_XX(iunif, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, du2, i0, i0, i1, i1, tmp, wrk1d, wrk2d, wrk3d)

  OPEN(20,file='partial.dat')
  DO i = 1,imax
     WRITE(20,'(7e)') x(i), dx(i,1), dx(i,2), u(i), due(i), du1(i), du2(i)
  ENDDO
  CLOSE(20)

  STOP
END PROGRAM VPARTIAL2
