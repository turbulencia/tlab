!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Stretching function s.t. the grid spacing dx/ds is given by a 
!# hyperbolic tangent function 
!#
!# dx/ds = f_0 + [(f_1-f_0)/2]*[ 1 + TANH[(s-s_1)/(2*delta_1)] ]
!#             + [(f_2-f_0)/2]*[ 1 + TANH[(s-s_2)/(2*delta_2)] ]
!#
!# Assumes just 1 segment
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE BLD_TANH(idir, x, imax, scalex)

  USE GRID_LOCAL

  IMPLICIT NONE

  TINTEGER idir, imax
  TREAL x(imax), scalex

! -----------------------------------------------------------------------
  TREAL s_1, f_1, delta_1
  TREAL s_2, f_2, delta_2
  TINTEGER i, iloc

! #######################################################################
! define parameters
  s_1    = iseg_vals(1,1,idir) ! transition point in uniform grid
  f_1    = iseg_vals(2,1,idir) ! ratio f_1/f_0
  delta_1= iseg_vals(3,1,idir)
  IF ( iseg_opts(2,1,idir) .EQ. 2 ) THEN ! mode 2
     s_2    = iseg_vals(4,1,idir) ! transition point in uniform grid
     f_2    = iseg_vals(5,1,idir) ! ratio f_2/f_0
     delta_2= iseg_vals(6,1,idir)
  ENDIF

! mirrowing case; first point in array is imax/2 
  IF ( idir_opts(3,idir) .EQ. 1 ) THEN; iloc = imax/2 ! mirrored case
  ELSE;                                 iloc = 1;     ENDIF

! create uniform reference grid s
  DO i=iloc,imax
     x(i) = M_REAL(i-iloc)/M_REAL(imax-iloc)*scalex
  ENDDO

! create grid x as a function of variable s 
  IF ( iseg_opts(2,1,idir) .EQ. 2 ) THEN
     DO i = iloc,imax
        x(i) = x(i) + (f_1-C_1_R)*delta_1*LOG(EXP((x(i)-s_1)/delta_1)+C_1_R) &
                    + (f_2-C_1_R)*delta_2*LOG(EXP((x(i)-s_2)/delta_2)+C_1_R)
     ENDDO
  ELSE
     DO i = iloc,imax
        x(i) = x(i) + (f_1-C_1_R)*delta_1*LOG(EXP((x(i)-s_1)/delta_1)+C_1_R)
     ENDDO
  ENDIF

! correct value of scale
  x(:) = x(:) - x(iloc)
  scalex = x(imax)

  RETURN
END SUBROUTINE BLD_TANH
