#include "types.h"

SUBROUTINE BLD_EXP(idir, x, imax, scalex)
  USE GRID_LOCAL
  IMPLICIT NONE

  TINTEGER, INTENT(IN   ) :: idir, imax
  TREAL,    INTENT(INOUT) :: x(imax), scalex

  ! -----------------------------------------------------------------------
  TREAL ds, s_old
  TINTEGER i, iloc, GRID_EXP, iseg
  EXTERNAL GRID_EXP

  ! quadrature variables
  TREAL ABSERR,EPSABS,EPSREL,RESULT
  TINTEGER IER,NEVAL
  TREAL WORK, params(9)
  TINTEGER IWORK,KEY,LAST,LENW,LIMIT
  DIMENSION IWORK(100),WORK(400)

  ! Convergence constants
  EPSABS = C_1EM10_R
  EPSREL = C_1EM10_R
  KEY = 6
  LIMIT = 100
  LENW = LIMIT*4

  ! #######################################################################
  iseg = 1                      ! Assumes just 1 segment

  iloc = 1
  IF ( g_build(idir)%mirrored ) iloc = imax/2 ! mirrored case; first point in array is imax/2

  ! create uniform reference grid s
  DO i=iloc,imax
    x(i) = M_REAL(i-iloc) /M_REAL(imax-iloc) *scalex
  ENDDO

  ! -----------------------------------------------------------------------
  ! define parameters; superposition of up to 3 modes, each with 3 parameters
  params(1:9) = g_build(idir)%vals(1:9,iseg)
  ds = scalex /M_REAL(imax-iloc)
  params(2::3) = params(2::3) /ds ! normalized stretching factor

  ! create mapping from grid s to grid x
  s_old   = C_0_R
  x(iloc) = C_0_R
  DO i=iloc+1,imax
    CALL QUADAG(GRID_EXP,params,s_old,x(i),EPSABS,&
        EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
    s_old  = x(i)
    x(i) = x(i-1) + RESULT
  ENDDO

  ! correct value of scale
  scalex = x(imax)

  RETURN
END SUBROUTINE BLD_EXP

! ###################################################################
TREAL FUNCTION GRID_EXP(y,p)
  TREAL, INTENT(IN) :: y, p(9)

  ! -----------------------------------------------------------------------
  TREAL st(3), df(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters

  ! #######################################################################
  ! define local variables for readability below; strides of 3
  st    = p(1::3)     ! transition point in uniform grid
  df    = p(2::3)     ! maximum stretching factor
  delta = p(3::3)

  GRID_EXP = C_1_R
  DO im = 1,3                           ! 3 modes at most
    IF ( ABS(delta(im)) > C_0_R ) THEN
      GRID_EXP = GRID_EXP *( EXP((y-st(im))/delta(im)) +C_1_R )**(df(im)*delta(im))
    END IF
  END DO
  !# Not sure why the previous version had
  !# GRID_EXP = GRID_EXP*( C_1_R + fds_3 /(COSH((y-s0_3)/delta_3))**C_2_R )

  RETURN
END FUNCTION GRID_EXP
