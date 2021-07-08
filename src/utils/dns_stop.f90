#include "types.h"

SUBROUTINE DNS_STOP(ic)

  USE DNS_CONSTANTS, ONLY : lfile, efile

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_time_min, ims_time_max, ims_time_trans, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  INTEGER, INTENT(IN) :: ic ! error code

  ! -------------------------------------------------------------------
  CHARACTER*256 line

  ! ###################################################################
  IF ( ic /= 0 ) THEN
    WRITE(line,*) ic
    line = 'Error code '//TRIM(ADJUSTL(line))//'.'
    CALL IO_WRITE_ASCII(efile,line)
  ENDIF

  CALL GETARG(0,line)
  WRITE(line,*) 'Finalizing program '//TRIM(ADJUSTL(line))
  IF ( ic == 0 ) THEN
    line = TRIM(ADJUSTL(line))//' normally.'
  ELSE
    line = TRIM(ADJUSTL(line))//' abnormally. Check '//TRIM(ADJUSTL(efile))
  ENDIF
  CALL IO_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
  ims_time_max = MPI_WTIME()
  WRITE(line,1000) ims_time_max-ims_time_min
  line = 'Time elapse ....................: '//TRIM(ADJUSTL(line))
  CALL IO_WRITE_ASCII(lfile, line)

#ifdef PROFILE_ON
  WRITE(line,1000) ims_time_trans
  line = 'Time in array transposition ....: '//TRIM(ADJUST(line))
  CALL IO_WRITE_ASCII(lfile, line)
#endif

1000 FORMAT(G_FORMAT_R)

#endif

  CALL IO_WRITE_ASCII(lfile, '########################################')

#ifdef USE_MPI
  CALL MPI_FINALIZE(ims_err)
#endif

  STOP
END SUBROUTINE DNS_STOP
