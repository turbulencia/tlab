PROGRAM VEFILTER
  
  IMPLICIT NONE

#include "types.h"
#include "integers.h"
  
  TINTEGER imode_fdm, imax, jmax, kmax, i, wk, i1bc, idummy, iunif
  PARAMETER(imax=256)
  TREAL scalex
  TREAL x(imax), dx(imax,2+4*3+4*3), u(imax), uf(imax), a(imax,5)
  TREAL wrk1d(imax,5), wrk2d(imax), wrk3d(imax)
  TREAL tmp(imax)
  
! ###################################################################
  scalex = C_2_R*C_PI_R
  jmax  = 1
  kmax  = 1
  imode_fdm = 6
  iunif = 1
  
  WRITE(*,*) 'Periodic (0) or nonperiodic (1) case ?'
  READ(*,*) i1bc
  
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
  
  ! TO BE REVIEWED
  ! CALL FDM_INITIALIZE(iunif, imode_fdm, imax, i1bc, scalex, x, dx, wrk1d)
  
! ###################################################################
! Define the function
  IF ( i1bc .EQ. 0 ) THEN
     WRITE(*,*) 'Wavenumber ?'
     READ(*,*) wk
     DO i = 1,imax
        u(i)  = SIN(C_2_R*C_PI_R/scalex*M_REAL(wk)*x(i) + C_PI_R*C_05_R)
        uf(i) = u(i)
     ENDDO
  ELSE
     OPEN(21,file='f.dat')
     DO i = 1,imax
        READ(21,*) u(i)
        uf(i) = u(i)
     ENDDO
     CLOSE(21)
  ENDIF

! ###################################################################
  ! CALL FILT4E_INI(imax, i1bc, scalex, x, a)
  ! CALL OPR_FILTER_X_OLD(i3, imax, jmax, kmax, i1bc, i0, i0, uf, a, wrk3d, wrk3d)
!  CALL  OPR_FILTER(i3, imax, jmax, kmax, kmax, i1bc, i1bc, i1bc, i1, i0, i0, i1, uf, a, a, a, wrk3d)

  OPEN(20,file='filter.dat')
  DO i = 1,imax
     WRITE(20,*) x(i), u(i), uf(i)
  ENDDO
  CLOSE(20)
  
  STOP
END PROGRAM VEFILTER
      
