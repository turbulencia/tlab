#include "types.h"

SUBROUTINE GRID_MIRROR(iflag, imax, x, scale)

  IMPLICIT NONE

  TINTEGER imax, iflag
  TREAL x(imax), scale

  TINTEGER i
  TREAL offset

! ####################################
! # Offset for even number of points #
! ####################################

  offset = (x(imax/2+1)-x(imax/2)) / C_2_R
  DO i = imax/2,imax
     x(i) = x(i) - offset
  ENDDO

! #############
! # Mirroring #
! #############

  DO i = 1,imax/2-1
     x(i) = - x(imax+1-i)
  ENDDO

! ######################
! # Global translation #
! ######################

  offset = x(1)
  DO i = 1,imax
     x(i) = x(i) - offset
  ENDDO

! #################
! # Final scaling #
! #################

  IF ( iflag .EQ. 1 ) THEN
     scale = C_2_R*scale

     DO i = 1,imax
        x(i) = x(i)/x(imax)*scale
     ENDDO

  ELSE IF ( iflag .EQ. 2 ) THEN
     scale = x(imax)-x(1)

  ENDIF

  RETURN
END SUBROUTINE GRID_MIRROR

