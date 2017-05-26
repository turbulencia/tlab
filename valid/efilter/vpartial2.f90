#include "types.h"

PROGRAM VPARTIAL2

  USE DNS_TYPES, ONLY : grid_dt
  IMPLICIT NONE

#include "integers.h"

  TYPE(grid_dt) :: g
  TINTEGER imax,jmax,kmax, i, wk, idummy
  PARAMETER(imax=128)
  TREAL x(imax,3+4*3+4*3), u(imax), du1(imax), du2(imax), due(imax)
  TREAL wrk1d(imax,5), wrk2d(imax), wrk3d(imax), bcs(2,2)
  TREAL tmp(imax)

! ###################################################################
  g%size     = imax
  g%scale    = C_2_R*C_PI_R
  g%mode_fdm = FDM_COM6_JACOBIAN
  g%uniform  = .TRUE.
  jmax = 1
  kmax = 1

  WRITE(*,*) 'Periodic (.TRUE. or .FALSE.)?'
  READ(*,*) g%periodic
  WRITE(*,*) 'Wavenumber ?'
  READ(*,*) wk

! CHANGE TO UPDATE NEW GRID_DT
  IF ( g%periodic ) THEN
     DO i = 1,imax
        g%nodes(i) = M_REAL(i-1)/M_REAL(imax)*g%scale
     ENDDO
  ELSE
     OPEN(21,file='y.dat')
     DO i = 1,imax
!        g%nodes(i) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
        READ(21,*) idummy, g%nodes(i)
     ENDDO
     CLOSE(21)
     g%scale = g%nodes(imax)-g%nodes(1)
  ENDIF

  CALL FDM_INITIALIZE(x, g, wrk1d)

! ###################################################################
! Define the function
  DO i = 1,imax
     u(i) = SIN(C_2_R*C_PI_R/g%scale*M_REAL(wk)*g%nodes(i))
     due(i) = -(C_2_R*C_PI_R/g%scale*M_REAL(wk))**2*u(i)
!     u(i) = EXP(-(g%nodes(i)-C_PI_R)**2/(C_PI_R**2/64.))
!     due(i) = -C_2_R*(g%nodes(i)-C_PI_R)/(C_PI_R**2/64.)*u(i)
  ENDDO

! ###################################################################
  bcs(:,1) = 0
  bcs(:,2) = 1
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs(1,1), g, u,   tmp, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs(1,2), g, tmp, du1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs,      g, u,   du2, tmp,   wrk2d,wrk3d)

  OPEN(20,file='partial.dat')
  DO i = 1,imax
     WRITE(20,'(7e)') g%nodes(i), g%jac(i,1), g%jac(i,2), u(i), due(i), du1(i), du2(i)
  ENDDO
  CLOSE(20)

  STOP
END PROGRAM VPARTIAL2
