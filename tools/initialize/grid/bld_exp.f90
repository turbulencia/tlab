#include "types.h"

SUBROUTINE BLD_EXP(idir, x, imax, scalex)
  
  USE GRID_LOCAL

  IMPLICIT NONE

  TINTEGER idir, imax
  TREAL x(imax), scalex

  TREAL ds, s_old
  TINTEGER i, iloc, GRID_EXP
  EXTERNAL GRID_EXP

! quadrature variables
  TREAL ABSERR,EPSABS,EPSREL,RESULT
  TINTEGER IER,NEVAL
  TREAL WORK, params(10)
  TINTEGER IWORK,KEY,LAST,LENW,LIMIT
  DIMENSION IWORK(100),WORK(400)

! mirrowing case; first point in array is imax/2 
  IF ( idir_opts(3,idir) .EQ. 1 ) THEN; iloc = imax/2 ! mirrored case
  ELSE;                                 iloc = 1;     ENDIF

! create uniform reference grid s
  DO i=iloc,imax
     x(i) = M_REAL(i-iloc)/M_REAL(imax-iloc)*scalex
  ENDDO

  ds = scalex/M_REAL(imax-iloc)

! Convergence constants
  EPSABS = C_1EM10_R
  EPSREL = C_1EM10_R
  KEY = 6
  LIMIT = 100
  LENW = LIMIT*4

! define parameters first segment
  params(1) = iseg_vals(2,1,idir)/ds ! f_1
  params(2) = iseg_vals(1,1,idir)    ! transition point in uniform grid s_1
  params(3) = iseg_vals(3,1,idir)    ! delta_1
  IF ( iseg_opts(2,1,idir) .EQ. 2 ) THEN
     params(4) = iseg_vals(5,1,idir)/ds ! f_2
     params(5) = iseg_vals(4,1,idir)    ! transition point in uniform grid s_2
     params(6) = iseg_vals(6,1,idir)    ! delta_2
  ELSE
     params(4) = C_0_R
  ENDIF

! create grid x as a function of variable s
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
FUNCTION GRID_EXP(y,p)
  
  TREAL y, p(*)
  TREAL GRID_EXP

  TREAL fds_1, s0_1, delta_1
  TREAL fds_2, s0_2, delta_2

! first segment
  fds_1   = p(1)
  s0_1    = p(2)
  delta_1 = p(3)

  GRID_EXP = (EXP((y-s0_1)/delta_1)+C_1_R)**(fds_1*delta_1)

! second segment
  fds_2   = p(4)
  s0_2    = p(5)
  delta_2 = p(6)

  IF ( fds_2 .GT. C_0_R ) THEN
     GRID_EXP = GRID_EXP*(EXP((y-s0_2)/delta_2)+C_1_R)**(fds_2*delta_2)
  ENDIF

  RETURN
END FUNCTION GRID_EXP
        
