#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 1996/07/26 - S. Stanley
!#              Created
!# 2010/03/28 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################

!########################################################################
! Read routine
!########################################################################
SUBROUTINE IO_READ_GRID(name, imax,jmax,kmax, scalex,scaley,scalez, x,y,z, area)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_PROCS

  IMPLICIT NONE

  CHARACTER*(*) name
  TINTEGER imax, jmax, kmax
  TREAL scalex, scaley, scalez
  TREAL x(imax), y(jmax), z(kmax)
  TREAL, OPTIONAL :: area

! -----------------------------------------------------------------------
  TINTEGER imaxdum, jmaxdum, kmaxdum
  TREAL scale(3)
  CHARACTER*(32) line

! #######################################################################
  OPEN(50,file=name, status='old',form='unformatted')
  REWIND(50)

! -----------------------------------------------------------------------
  READ(50) imaxdum, jmaxdum, kmaxdum
  READ(50) scale
  scalex = scale(1)
  scaley = scale(2)
  scalez = scale(3)

! -----------------------------------------------------------------------
  IF (imaxdum .NE. imax .OR. jmaxdum .NE. jmax .OR. kmaxdum .NE. kmax) THEN
     CLOSE(50)
     WRITE(line,100) imaxdum,jmaxdum,kmaxdum
     CALL TLAB_WRITE_ASCII(efile, 'IO_READ_GRID. Dimensions ('//TRIM(line)//') unmatched.')
     CALL TLAB_STOP(DNS_ERROR_DIMGRID)
  ENDIF

! -----------------------------------------------------------------------
  READ(50) x
  READ(50) y
  READ(50) z
  CLOSE(50)

  IF ( PRESENT(area) ) THEN
    area = scalex
    IF ( kmax > 1 ) area = area *scalez ! 3D case
  ENDIF

  RETURN

100 FORMAT(I5,',',I5,',',I5)

END SUBROUTINE IO_READ_GRID

!########################################################################
! Write routine
!########################################################################
SUBROUTINE IO_WRITE_GRID(name, imax,jmax,kmax, scalex,scaley,scalez, x,y,z)

  IMPLICIT NONE

  CHARACTER*(*) name
  TINTEGER imax, jmax, kmax
  TREAL scalex, scaley, scalez
  TREAL x(imax), y(jmax), z(kmax)

! -----------------------------------------------------------------------
  TREAL scale(3)

!########################################################################
  OPEN(unit=51,file=name,form='unformatted', status='unknown')

! -----------------------------------------------------------------------
  scale(1) = scalex
  scale(2) = scaley
  scale(3) = scalez
  WRITE(51) imax, jmax, kmax
  WRITE(51) scale

! -----------------------------------------------------------------------
  WRITE(51) x
  WRITE(51) y
  WRITE(51) z
  CLOSE(51)

  RETURN
END SUBROUTINE IO_WRITE_GRID
