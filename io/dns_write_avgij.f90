#include "types.h"

#define LOC_UNIT_ID i57
#define LOC_STATUS 'unknown'

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2001/08/01 - J.P. Mellado
!#              Created
!# 2008/01/09 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Note that the array mean1d is duplicated in each processor, and only
!# that from PE0 is written
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE DNS_WRITE_AVGIJ(name, mean1d)

  USE DNS_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro
#endif
#ifdef LES
  USE LES_GLOBAL
#endif

  IMPLICIT NONE

#include "integers.h"

  CHARACTER*(*) name
  TREAL mean1d(nstatavg,jmax,inb_mean_spatial)

! -------------------------------------------------------------------
  CHARACTER*128 :: line

! ###################################################################
  line = 'Writing field '//TRIM(ADJUSTL(name))//'...'

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
#include "dns_open_file.h"
     CALL WRT_STHD(LOC_UNIT_ID, i0, istat_maj_ver, istat_min_ver, &
          itime, rtime, istattimeorg, rstattimeorg,&
          nstatavg, jmax, inb_mean_spatial, nstatavg_points, statavg)
#ifdef LES
     IF ( iles .EQ. 1 ) THEN
        WRITE(LOC_UNIT_ID) ilesstat_maj_ver, ilesstat_min_ver, iarmavg_pts
     ENDIF
#endif
     WRITE(LOC_UNIT_ID) mean1d
     CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE DNS_WRITE_AVGIJ

! ###################################################################
! ###################################################################
SUBROUTINE WRT_STHD(unit, irec, major, minor, &
     iter, rtime, iterorg, rtimeorg,&
     nstatavg, jmax, nstat, nstatavg_points, statavg)

  IMPLICIT NONE

  TINTEGER unit, irec
  TINTEGER iter, iterorg
  TINTEGER nstatavg, jmax, nstat, nstatavg_points
  TREAL rtime, rtimeorg
  TINTEGER statavg(nstatavg)
  TINTEGER major, minor

  TREAL tmp(1)
  TINTEGER reclen

  reclen = SIZEOFINT*2
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) major
     WRITE(unit) minor
     WRITE(unit) reclen
  ELSE
     WRITE(unit) major, minor
  ENDIF

  tmp(1) = rtime
  reclen = SIZEOFINT+SIZEOFREAL 
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) iter
     WRITE(unit) rtime
     WRITE(unit) reclen
  ELSE
     WRITE(unit) iter, tmp
  ENDIF

  tmp(1) = rtimeorg
  reclen = SIZEOFINT+SIZEOFREAL 
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) iterorg
     WRITE(unit) rtimeorg
     WRITE(unit) reclen
  ELSE
     WRITE(unit) iterorg, tmp
  ENDIF

  reclen = 4*SIZEOFINT
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) nstatavg
     WRITE(unit) jmax
     WRITE(unit) nstat
     WRITE(unit) nstatavg_points
     WRITE(unit) reclen
  ELSE
     WRITE(unit) nstatavg, jmax, nstat, nstatavg_points
  ENDIF

  reclen = nstatavg*SIZEOFINT
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) statavg
     WRITE(unit) reclen
  ELSE
     WRITE(unit) statavg
  ENDIF

  RETURN
END SUBROUTINE WRT_STHD

