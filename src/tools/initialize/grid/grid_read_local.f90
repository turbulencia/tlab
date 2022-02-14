#include "types.h"
#include "dns_error.h"

SUBROUTINE GRID_READ_LOCAL (inifile, idir, scale, periodic)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_PROCS
  USE GRID_LOCAL

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)    :: inifile
  TINTEGER,         INTENT(IN)    :: idir
  TREAL,            INTENT(  OUT) :: scale
  LOGICAL,          INTENT(  OUT) :: periodic

! -------------------------------------------------------------------
  CHARACTER*64 sRes, str, bakfile
  CHARACTER*9 title
  TINTEGER iseg, idummy

! #######################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

! Defining direction
  IF      ( idir == 1 ) THEN; title = 'IniGridOx'
  ELSE IF ( idir == 2 ) THEN; title = 'IniGridOy'
  ELSE IF ( idir == 3 ) THEN; title = 'IniGridOz'
  ENDIF

  CALL TLAB_WRITE_ASCII(bakfile, '['//title//']')

! -------------------------------------------------------------------
! General options
! -------------------------------------------------------------------
  CALL TLAB_WRITE_ASCII(bakfile, 'segments=<number of segments>')
  CALL TLAB_WRITE_ASCII(bakfile, 'periodic=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, 'mirrored=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, 'fixed_scale=<value>')

  CALL SCANINIINT(bakfile, inifile, title, 'segments', '1', g_build(idir)%nseg)

  periodic = .FALSE.
  CALL SCANINICHAR(bakfile, inifile, title, 'periodic', 'no', sRes)
  IF ( TRIM(ADJUSTL(sRes)) == 'yes' ) periodic = .TRUE.

  g_build(idir)%mirrored = .FALSE.
  CALL SCANINICHAR(bakfile, inifile, title, 'mirrored', 'no', sRes)
  IF ( TRIM(ADJUSTL(sRes)) == 'yes' ) g_build(idir)%mirrored = .TRUE.

  CALL SCANINIREAL(bakfile, inifile, title, 'fixed_scale', '-1.0', g_build(idir)%fixed_scale )

  IF ( periodic .AND. g_build(idir)%mirrored ) THEN
     CALL TLAB_WRITE_ASCII(efile, 'GRID_READ_LOCAL. Periodicity with mirroring is not supported.')
     CALL TLAB_STOP(DNS_ERROR_GRID_SCALE)
  ENDIF

! -------------------------------------------------------------------
! Loop over the segments
! -------------------------------------------------------------------
  DO iseg = 1,g_build(idir)%nseg
     WRITE(str,*) iseg

     CALL TLAB_WRITE_ASCII(bakfile, 'Segment number '//TRIM(ADJUSTL(str)) )
     CALL TLAB_WRITE_ASCII(bakfile, 'scales_'//TRIM(ADJUSTL(str))//'=<physical end of the segment>')
     CALL TLAB_WRITE_ASCII(bakfile, 'points_'//TRIM(ADJUSTL(str))//'=<points in the segment>')
     CALL TLAB_WRITE_ASCII(bakfile, 'opts_'//TRIM(ADJUSTL(str))//'=<option>')
     CALL TLAB_WRITE_ASCII(bakfile, 'vals_'//TRIM(ADJUSTL(str))//'=<values>')


     CALL SCANINIINT(bakfile,  inifile, title, 'points_'//TRIM(ADJUSTL(str)),   '1', g_build(idir)%size(iseg) )
     CALL SCANINIREAL(bakfile, inifile, title, 'scales_'//TRIM(ADJUSTL(str)),'-1.0', g_build(idir)%end(iseg)  )

     g_build(idir)%opts(:,iseg) = 0
     CALL SCANINICHAR(bakfile, inifile, title, 'opts_'//TRIM(ADJUSTL(str)), '1', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) == 'uniform' ) THEN; g_build(idir)%opts(1,iseg) = GTYPE_UNIFORM
     ELSE IF ( TRIM(ADJUSTL(sRes)) == 'tanh'    ) THEN; g_build(idir)%opts(1,iseg) = GTYPE_TANH
     ELSE IF ( TRIM(ADJUSTL(sRes)) == 'exp'     ) THEN; g_build(idir)%opts(1,iseg) = GTYPE_EXP
     ELSE
       idummy = MAX_PARAMES
       CALL LIST_INTEGER(sRes, idummy, g_build(idir)%opts(1,iseg))
     END IF

     g_build(idir)%vals(:,iseg) = 0
     CALL SCANINICHAR(bakfile, inifile, title, 'vals_'//TRIM(ADJUSTL(str)), '1.0', sRes)
     idummy = MAX_PARAMES
     CALL LIST_REAL(sRes, idummy, g_build(idir)%vals(1,iseg))

  ENDDO

! -------------------------------------------------------------------
! Control
! -------------------------------------------------------------------
  scale = g_build(idir)%end(g_build(idir)%nseg)
  IF ( scale <= C_0_R ) THEN
     CALL TLAB_WRITE_ASCII(efile, 'GRID_READ_LOCAL. Scales undefined.')
     CALL TLAB_STOP(DNS_ERROR_GRID_SCALE)

  ENDIF

  RETURN
END SUBROUTINE GRID_READ_LOCAL
