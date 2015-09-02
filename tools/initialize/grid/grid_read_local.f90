#include "types.h"
#include "dns_error.h"

SUBROUTINE GRID_READ_LOCAL (inifile,idir)

  USE DNS_CONSTANTS, ONLY : efile
  USE GRID_LOCAL

  IMPLICIT NONE

  CHARACTER*(*) inifile
  TINTEGER idir

! -------------------------------------------------------------------
  CHARACTER*64 sRes, str, bakfile
  CHARACTER*9 title
  TINTEGER iseg, idummy

! #######################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

! Defining direction
  IF      (idir.eq.1) THEN; title = 'IniGridOx'
  ELSE IF (idir.eq.2) THEN; title = 'IniGridOy'
  ELSE IF (idir.eq.3) THEN; title = 'IniGridOz'; ENDIF

  CALL IO_WRITE_ASCII(bakfile, '['//title//']')

! -------------------------------------------------------------------
! General options
! -------------------------------------------------------------------
  CALL IO_WRITE_ASCII(bakfile, 'segments=<number of segments>')
  CALL IO_WRITE_ASCII(bakfile, 'periodic=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, 'mirrored=<yes/no>')

  CALL SCANINIINT(bakfile, inifile, title, 'segments', '1', idir_opts(1,idir))

  CALL SCANINICHAR(bakfile, inifile, title, 'periodic', 'no', sRes) 
  IF (TRIM(ADJUSTL(sRes)) .eq. 'yes') THEN; idir_opts(2,idir) = 1
  ELSE;                                     idir_opts(2,idir) = 0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, title, 'mirrored', 'no', sRes)
  IF (TRIM(ADJUSTL(sRes)) .eq. 'yes') THEN; idir_opts(3,idir) = 1
  ELSE;                                     idir_opts(3,idir) = 0; ENDIF

  IF (idir_opts(2,idir).eq.1 .and. idir_opts(3,idir).eq.1) THEN
     CALL IO_WRITE_ASCII(efile, 'GRID_READ_LOCAL. Periodicity with mirroring is not supported.')
     CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
  ENDIF

! -------------------------------------------------------------------
! Loop over the segments
! -------------------------------------------------------------------
  DO iseg = 1,idir_opts(1,idir)
     WRITE(str,*) iseg

     CALL IO_WRITE_ASCII(bakfile, 'Segment number '//TRIM(ADJUSTL(str)) )
     CALL IO_WRITE_ASCII(bakfile, 'scales_'//TRIM(ADJUSTL(str))//'=<physical end of the segment>')
     CALL IO_WRITE_ASCII(bakfile, 'points_'//TRIM(ADJUSTL(str))//'=<points in the segment>')
     CALL IO_WRITE_ASCII(bakfile, 'opts_'//TRIM(ADJUSTL(str))//'=<option>')
     CALL IO_WRITE_ASCII(bakfile, 'vals_'//TRIM(ADJUSTL(str))//'=<values>')


     CALL SCANINIINT(bakfile,  inifile, title, 'points_'//TRIM(ADJUSTL(str)),'1', isegdim(iseg,idir))
     iseg_opts(:,iseg,idir) = 0
     CALL SCANINICHAR(bakfile, inifile, title, 'opts_'//TRIM(ADJUSTL(str)), '1', sRes)
     idummy = MAX_PARAMES
     CALL LIST_INTEGER(sRes, idummy, iseg_opts(1,iseg,idir))

     CALL SCANINIREAL(bakfile, inifile, title, 'scales_'//TRIM(ADJUSTL(str)),'-1.0', isegend(iseg,idir))
     iseg_vals(:,iseg,idir) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, title, 'vals_'//TRIM(ADJUSTL(str)), '1.0', sRes)
     idummy = MAX_PARAMES
     CALL LIST_REAL(sRes, idummy, iseg_vals(1,iseg,idir))

  ENDDO

! -------------------------------------------------------------------
! Control
! -------------------------------------------------------------------
  scale(idir) = isegend(idir_opts(1,idir),idir)
  IF ( scale(idir) .LE. C_0_R ) THEN
     CALL IO_WRITE_ASCII(efile, 'GRID_READ_LOCAL. Scales undefined.')
     CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
     
  ENDIF

  RETURN
END SUBROUTINE GRID_READ_LOCAL

