#include "types.h"
#include "dns_error.h"
  
!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_READ_CHEMKIN(name)

  USE THERMO_GLOBAL
  USE DNS_CONSTANTS, ONLY : efile, lfile

  IMPLICIT NONE

#include "integers.h"

  CHARACTER*(*) name

! -----------------------------------------------------------------------
  TREAL T1, T2, T3
  TINTEGER i, il, is
  LOGICAL frun
  CHARACTER*15 token
  CHARACTER*225 wline
  CHARACTER*80 line, line1, line2, line3
  TINTEGER THERMO_FLAG(MAX_NSP)

! #######################################################################
! Initialize thermodynamic data structure
  DO is=1, NSP
     THERMO_FLAG(is) = 0
  ENDDO

! Read Thermodynamic file
  OPEN(i23,file=name,status='old')

  REWIND(i23)

! Read Header
  READ(i23,*) line
  CALL IO_WRITE_ASCII(lfile, line)

  IF ( TRIM(ADJUSTL(line)) .NE. 'THERMO' ) THEN
     CALL IO_WRITE_ASCII(efile, 'THERMO_READ_CHEMKIN. Thermodynamic file format error')
     CALL DNS_STOP(DNS_ERROR_THERMOFORMAT)
  ENDIF

! Read Temperature ranges
  READ(i23, *) T1, T2, T3
  WRITE(wline,*) T1, T2, T3
  CALL IO_WRITE_ASCII(lfile, wline(1:80))

! Remove comments
  frun = .TRUE.
  DO WHILE( frun )
     READ(i23,'(A80)',END=50) line
     IF ( line(1:1) .NE. '!' ) frun = .FALSE.
  ENDDO

  frun = .TRUE. 
  DO WHILE ( frun )
!    Check for end of file
     READ(line,'(A15)',END=50) token
     IF ( TRIM(ADJUSTL(token)) .EQ. 'END' ) THEN
        frun = .FALSE.
        GOTO 50
     ENDIF

!    Read all relevant information
     READ(i23,'(A80)',END=50) line1
     READ(i23,'(A80)',END=50) line2
     READ(i23,'(A80)',END=50) line3

!    Process lines
     DO is = 1,NSP
        IF ( TRIM(ADJUSTL(token)) .EQ. THERMO_SPNAME(is) ) THEN
           CALL IO_WRITE_ASCII(lfile, line)
           CALL IO_WRITE_ASCII(lfile, line1)
           CALL IO_WRITE_ASCII(lfile, line2)
           CALL IO_WRITE_ASCII(lfile, line3)

!          Required species found, process information
!          Get limit temperatures
           DO i=1, 225
              wline(i:i) = ' '
           ENDDO
           wline=line(46:75)
           READ(wline,*) (THERMO_TLIM(i,is),i=1,3)

!          Concatenate lines so read is simpler
           wline=line1(1:75)//line2(1:75)//line3(1:75)

           DO i=1,14
              il = (i-1)*15+1
              READ(wline(il:il+14),*) THERMO_AI(i,1,is)
           ENDDO

           THERMO_FLAG(is) = 1
        ENDIF
     ENDDO

!    Read next line
     READ(i23, '(A80)',END=50) line

  ENDDO

50 CLOSE(i23)

  DO is = 1,NSP
     IF ( THERMO_FLAG(is) .EQ. 0 ) THEN
        CALL IO_WRITE_ASCII(efile, 'THERMO_READ_CHEMKIN. Not all thermodynamic data contained in thermo file')
        CALL DNS_STOP(DNS_ERROR_THERMOCONT)
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE THERMO_READ_CHEMKIN
