#include "types.h"
#include "dns_error.h"

#define LOC_UNIT_ID i56
#define LOC_STATUS 'old'

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
!# that from PE0 is read
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE DNS_READ_AVGIJ(name,mean1d)
  
  USE DNS_GLOBAL, ONLY : istat_maj_ver, istat_min_ver
  USE DNS_GLOBAL, ONLY : istattimeorg, rstattimeorg, nstatavg_points, nstatavg, statavg
  USE DNS_GLOBAL, ONLY : itime, rtime, jmax, inb_mean_spatial
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
  TREAL, DIMENSION(nstatavg*jmax*inb_mean_spatial) :: mean1d
  
! -------------------------------------------------------------------
  CHARACTER*128 :: line
  TINTEGER nstat
  LOGICAL lfilexist

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
        nstat = inb_mean_spatial
        CALL RD_STHD(LOC_UNIT_ID, i0, istat_maj_ver, istat_min_ver, &
             itime, rtime, istattimeorg, rstattimeorg, &
             nstatavg, jmax, nstat, nstatavg_points, statavg)
#ifdef LES
        IF ( iles .EQ. 1 ) THEN
           READ(LOC_UNIT_ID) ilesstat_maj_ver, ilesstat_min_ver, iarmavg_pts
        ENDIF
#endif
        READ(LOC_UNIT_ID) mean1d(1:nstatavg*jmax*nstat)
        CLOSE(LOC_UNIT_ID)

! -------------------------------------------------------------------
! Initialize data
! -------------------------------------------------------------------
     ELSE
        nstatavg_points = 0
        istattimeorg = itime
        rstattimeorg = rtime
#ifdef LES
        IF ( iles .EQ. 1 ) THEN
           ilesstat_maj_ver = 1
           ilesstat_min_ver = 2
           iarmavg_pts = 0
        ENDIF
#endif
        mean1d = C_0_R
        CALL IO_WRITE_ASCII(lfile,'Statistics have been initialized.')
     ENDIF
     
#ifdef USE_MPI
  ENDIF
#endif
  
  RETURN
END SUBROUTINE DNS_READ_AVGIJ

! ###################################################################
! ###################################################################
SUBROUTINE RD_STHD(unit, irec, major, minor, iter, rtime, iterorg, rtimeorg,&
     nstatavg, jmax, nstat, nstatavg_points, statavg)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER unit, irec
  TINTEGER iter, iterorg
  TINTEGER nstatavg, jmax, nstat, nstatavg_points
  TREAL rtime, rtimeorg
  TINTEGER statavg(nstatavg)

  TINTEGER major, minor
  TINTEGER iterdum
  TINTEGER nstatavgdum, jmaxdum, nstatdum
  TREAL tmp(1)
  TINTEGER reclen

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) major
     READ(unit) minor
     READ(unit) reclen
  ELSE
     READ(unit) major, minor
  ENDIF

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

