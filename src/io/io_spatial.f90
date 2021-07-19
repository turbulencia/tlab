#include "types.h"
#include "dns_error.h"
#include "avgij_map.h"

#define LOC_UNIT_ID i57
#define LOC_STATUS 'unknown'

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
SUBROUTINE IO_WRITE_AVG_SPATIAL(name, mean_flow, mean_scal)

  USE DNS_GLOBAL, ONLY : istattimeorg, rstattimeorg, nstatavg_points, nstatavg, statavg
  USE DNS_GLOBAL, ONLY : itime, rtime, jmax, inb_scal
  USE DNS_CONSTANTS, ONLY : lfile

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro
#endif

  IMPLICIT NONE

#include "integers.h"

  CHARACTER*(*) name
  TREAL mean_flow(nstatavg,jmax,MA_MOMENTUM_SIZE)
  TREAL mean_scal(nstatavg,jmax,MS_SCALAR_SIZE,inb_scal)

! -------------------------------------------------------------------
  CHARACTER*128 :: line
  TINTEGER nstat

! ###################################################################
  line = 'Writing field '//TRIM(ADJUSTL(name))//'...'

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
#include "dns_open_file.h"
     nstat = MA_MOMENTUM_SIZE +MS_SCALAR_SIZE*inb_scal
     CALL WRT_STHD(LOC_UNIT_ID, i0, &
          itime, rtime, istattimeorg, rstattimeorg,&
          nstatavg, jmax, nstat, nstatavg_points, statavg)
     WRITE(LOC_UNIT_ID) mean_flow
     WRITE(LOC_UNIT_ID) mean_scal
     CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_WRITE_AVG_SPATIAL

! ###################################################################
! ###################################################################
SUBROUTINE WRT_STHD(unit, irec, &
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

#undef LOC_STATUS

!########################################################################
!########################################################################
#define LOC_STATUS 'old'

SUBROUTINE IO_READ_AVG_SPATIAL(name,mean_flow,mean_scal)

  USE DNS_GLOBAL, ONLY : istattimeorg, rstattimeorg, nstatavg_points, nstatavg, statavg
  USE DNS_GLOBAL, ONLY : itime, rtime, jmax, inb_scal
  USE DNS_CONSTANTS, ONLY : lfile

#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) name
  TREAL, DIMENSION(nstatavg*jmax*MA_MOMENTUM_SIZE)        :: mean_flow
  TREAL, DIMENSION(nstatavg*jmax*MS_SCALAR_SIZE*inb_scal) :: mean_scal

! -------------------------------------------------------------------
  CHARACTER*128 :: line
  LOGICAL lfilexist
  TINTEGER nstat

! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     lfilexist = .FALSE.
     INQUIRE(file=name,EXIST=lfilexist)

! -------------------------------------------------------------------
! Read data
! -------------------------------------------------------------------
     IF ( lfilexist ) THEN
        line = 'Reading field '//TRIM(ADJUSTL(name))//'...'
        CALL IO_WRITE_ASCII(lfile,line)

#include "dns_open_file.h"
        REWIND(LOC_UNIT_ID)
        nstat = MA_MOMENTUM_SIZE +MS_SCALAR_SIZE*inb_scal
        CALL RD_STHD(LOC_UNIT_ID, i0, &
             itime, rtime, istattimeorg, rstattimeorg, &
             nstatavg, jmax, nstat, nstatavg_points, statavg)
#ifdef LES
        IF ( iles .EQ. 1 ) THEN
           READ(LOC_UNIT_ID) ilesstat_maj_ver, ilesstat_min_ver, iarmavg_pts
        ENDIF
#endif
        READ(LOC_UNIT_ID) mean_flow
        READ(LOC_UNIT_ID) mean_scal
        CLOSE(LOC_UNIT_ID)

! -------------------------------------------------------------------
! Initialize data
! -------------------------------------------------------------------
     ELSE
        nstatavg_points = 0
        istattimeorg = itime
        rstattimeorg = rtime
        mean_flow = C_0_R
        mean_scal = C_0_R
        CALL IO_WRITE_ASCII(lfile,'Statistics have been initialized.')
     ENDIF

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_READ_AVG_SPATIAL

! ###################################################################
! ###################################################################
SUBROUTINE RD_STHD(unit, irec, iter, rtime, iterorg, rtimeorg,&
     nstatavg, jmax, nstat, nstatavg_points, statavg)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER unit, irec
  TINTEGER iter, iterorg
  TINTEGER nstatavg, jmax, nstat, nstatavg_points
  TREAL rtime, rtimeorg
  TINTEGER statavg(nstatavg)

  TINTEGER iterdum
  TINTEGER nstatavgdum, jmaxdum, nstatdum
  TREAL tmp(1)
  TINTEGER reclen

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) iterdum
     READ(unit) rtime
     READ(unit) reclen
  ELSE
     READ(unit) iterdum, tmp
     rtime = tmp(1)
  ENDIF

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) iterorg
     READ(unit) rtimeorg
     READ(unit) reclen
  ELSE
     READ(unit) iterorg, tmp
     rtimeorg = tmp(1)
  ENDIF

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) nstatavgdum
     READ(unit) jmaxdum
     READ(unit) nstatdum
     READ(unit) nstatavg_points
     READ(unit) reclen
  ELSE
     READ(unit) nstatavgdum, jmaxdum, nstatdum, nstatavg_points
  ENDIF

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) statavg
     READ(unit) reclen
  ELSE
     READ(unit) statavg
  ENDIF

! #####################
! # Checking
! #####################

  IF (iterdum .NE. iter) THEN
     CALL IO_WRITE_ASCII(efile,'Stat file error (iter mismatch).')
     CALL DNS_STOP(DNS_ERROR_STFILE)
  ENDIF

  IF (jmaxdum .NE. jmax) THEN
     CALL IO_WRITE_ASCII(efile,'Stat file error (jmax mismatch).')
     CALL DNS_STOP(DNS_ERROR_STFILE)
  ENDIF

  IF (nstatavgdum .NE. nstatavg) THEN
     CALL IO_WRITE_ASCII(efile,'Stat file error (nstatavg mismatch).')
     CALL DNS_STOP(DNS_ERROR_STFILE)
  ENDIF

  IF (nstatdum .NE. nstat) THEN
     CALL IO_WRITE_ASCII(efile,'Stat file error (nstat mismatch).')
     CALL DNS_STOP(DNS_ERROR_STFILE)
  ELSE
     nstat = nstatdum
  ENDIF

  RETURN
END SUBROUTINE RD_STHD
