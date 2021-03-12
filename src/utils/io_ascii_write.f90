#include "types.h"

! #######################################################################
! Write ASCII data. Outside of ASCII library because of dependencies issues
! #######################################################################
SUBROUTINE IO_WRITE_ASCII(file, line)
  
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL, ONLY : imode_verbosity
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro
#endif  
  IMPLICIT NONE

  CHARACTER*(*), INTENT(IN) :: file, line

! -----------------------------------------------------------------------
  CHARACTER*10 clock(2)

! #######################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

     IF ( imode_verbosity .GT. 0 ) THEN

        OPEN(UNIT=22, FILE=file, STATUS='unknown',POSITION='APPEND')
        IF      ( imode_verbosity .EQ. 1 ) THEN
           WRITE(22,'(a)') TRIM(ADJUSTL(line))
        ELSE IF ( imode_verbosity .EQ. 2 ) THEN
           CALL DATE_AND_TIME(clock(1),clock(2))
           WRITE(22,'(a)') '['//TRIM(ADJUSTR(clock(2)))//'] '//TRIM(ADJUSTL(line))
        ENDIF
        CLOSE(22)

     ENDIF

#ifndef PARALLEL
     IF ( file .EQ. efile ) THEN
        WRITE(*,*) TRIM(ADJUSTL(line))
     ENDIF
#endif

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_WRITE_ASCII

! #######################################################################
! #######################################################################
SUBROUTINE IO_WRITE_ASCII_ALL(file, line)
  
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : imode_verbosity
  IMPLICIT NONE

  CHARACTER*(*), INTENT(IN) :: file, line

! -----------------------------------------------------------------------
  CHARACTER*10 clock(2)

! #######################################################################
  IF ( imode_verbosity .GT. 0 ) THEN
     
     OPEN(UNIT=22, FILE=file, STATUS='unknown',POSITION='APPEND')
     IF      ( imode_verbosity .EQ. 1 ) THEN
        WRITE(22,'(a)') TRIM(ADJUSTL(line))
     ELSE IF ( imode_verbosity .EQ. 2 ) THEN
        CALL DATE_AND_TIME(clock(1),clock(2))
        WRITE(22,'(a)') '['//TRIM(ADJUSTR(clock(2)))//'] '//TRIM(ADJUSTL(line))
     ENDIF
     CLOSE(22)

  ENDIF

#ifndef PARALLEL
  IF ( file .EQ. efile ) THEN
     WRITE(*,*) TRIM(ADJUSTL(line))
  ENDIF
#endif
  
  RETURN
END SUBROUTINE IO_WRITE_ASCII_ALL
