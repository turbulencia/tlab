#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2001/01/01 - C. Pantano
!#              Created
!# 2011/08/01 - J.P. Mellado
!#              Case insensitive
!#
!########################################################################
!# DESCRIPTION
!#
!# Scan ASCII file for inputs
!# Strings are always forced into lower-case
!#
!########################################################################

! #######################################################################
! Scan file for an integer value
! #######################################################################
SUBROUTINE SCANINIINT(ofile, ifile, title, name, default, value)
  IMPLICIT NONE
  
  CHARACTER*(*), INTENT(IN)  :: ofile, ifile, title, name, default
  TINTEGER,      INTENT(OUT) :: value

  CHARACTER*(128) StrValue

  CALL IO_READ_ASCII(ifile, title, name, StrValue, default)
  READ(StrValue, *) value
  
  CALL IO_WRITE_ASCII(ofile, TRIM(ADJUSTL(name))//'='//TRIM(ADJUSTL(StrValue)) )
  
  RETURN
END SUBROUTINE SCANINIINT

! #######################################################################
! Scan file for an integer value
! #######################################################################
SUBROUTINE SCANINILONGINT(ofile, ifile, title, name, default, value)
  IMPLICIT NONE
  
  CHARACTER*(*), INTENT(IN)  :: ofile, ifile, title, name, default
  TLONGINTEGER,  INTENT(OUT) :: value

  CHARACTER*(128) StrValue

  CALL IO_READ_ASCII(ifile, title, name, StrValue, default)
  READ(StrValue, *) value
  
  CALL IO_WRITE_ASCII(ofile, TRIM(ADJUSTL(name))//'='//TRIM(ADJUSTL(StrValue)) )
  
  RETURN
END SUBROUTINE SCANINILONGINT

! #######################################################################
! Scan file for an real value
! #######################################################################
SUBROUTINE SCANINIREAL(ofile, ifile, title, name, default, value)
  IMPLICIT NONE

  CHARACTER*(*), INTENT(IN)  :: ofile, ifile, title, name, default 
  TREAL,         INTENT(OUT) :: value

  CHARACTER*(128) StrValue

  CALL IO_READ_ASCII(ifile, title, name, StrValue, default)
  READ(StrValue, *) value

  CALL IO_WRITE_ASCII(ofile, TRIM(ADJUSTL(name))//'='//TRIM(ADJUSTL(StrValue)) )
  
  RETURN
END SUBROUTINE SCANINIREAL

! #######################################################################
! Scan file for an char value
! #######################################################################
SUBROUTINE SCANINICHAR(ofile, ifile, title, name, default, value)
  IMPLICIT NONE
  
  CHARACTER*(*), INTENT(IN)  :: ofile, ifile, title, name, default
  CHARACTER*(*), INTENT(OUT) :: value
  
  CALL IO_READ_ASCII(ifile, title, name, value, default)
  
  CALL IO_WRITE_ASCII(ofile, TRIM(ADJUSTL(name))//'='//TRIM(ADJUSTL(value)) )
  
  RETURN
END SUBROUTINE SCANINICHAR

! #######################################################################
! Scan file for a string
! #######################################################################
SUBROUTINE IO_READ_ASCII(fname, title, name, value, default)

#ifdef USE_MPI  
  USE DNS_MPI, ONLY : ims_pro, ims_err
#endif
  IMPLICIT NONE

#ifdef USE_MPI 
#include "mpif.h" 
#endif
  CHARACTER*(*), INTENT(IN)  :: fname, title, name, default
  CHARACTER*(*), INTENT(OUT) :: value

! -----------------------------------------------------------------------
  CHARACTER*512 line
  CHARACTER*128 str, tag1, tag2, tag3
  TINTEGER equal, code, n

! #######################################################################
! make case-insensitive: title, name and default tags.
  tag1=title; tag2=name; tag3=default
  DO n = 1,LEN(tag1)
     code = IACHAR(tag1(n:n)); IF ( code .GE. 65 .AND. code .LE. 90 ) tag1(n:n) = ACHAR(code+32)
  ENDDO
  DO n = 1,LEN(tag2)
     code = IACHAR(tag2(n:n)); IF ( code .GE. 65 .AND. code .LE. 90 ) tag2(n:n) = ACHAR(code+32)
  ENDDO
  DO n = 1,LEN(tag3)
     code = IACHAR(tag3(n:n)); IF ( code .GE. 65 .AND. code .LE. 90 ) tag3(n:n) = ACHAR(code+32)
  ENDDO

  value = tag3 !default

! -----------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN 
#endif 

  OPEN(unit=45,file=fname,status='OLD')

20 CONTINUE
  READ(45,'(A512)',END=50) line; line=TRIM(ADJUSTL(line))
  DO n = 1,LEN(line)
     code = IACHAR(line(n:n)); IF ( code .GE. 65 .AND. code .LE. 90 ) line(n:n) = ACHAR(code+32)
  ENDDO
  
  IF ( TRIM(ADJUSTL(line)) .EQ. '['//TRIM(ADJUSTL(tag1))//']' ) THEN ! Scan within block

30   CONTINUE
     READ(45,'(A512)',END=50) line; line=TRIM(ADJUSTL(line))
     IF ( line(1:1) .EQ. '[' ) GOTO 20 ! Scape sequence
     IF ( line(1:1) .EQ. '#' ) GOTO 30 ! Scape sequence
     equal = INDEX(line,'='); str = TRIM(ADJUSTL(line(1:equal-1)))
     DO n = 1,LEN(str)
        code = IACHAR(str(n:n)); IF ( code .GE. 65 .AND. code .LE. 90 ) str(n:n) = ACHAR(code+32)
     ENDDO

     IF ( str .EQ. tag2 ) THEN
        value = TRIM(ADJUSTL(line(equal+1:)))
        DO n = 1,LEN(value)
           code = IACHAR(value(n:n)); IF ( code .GE. 65 .AND. code .LE. 90 ) value(n:n) = ACHAR(code+32)
        ENDDO
        GOTO 50
     ENDIF
     GOTO 30

  ENDIF
  GOTO 20

50 CLOSE(unit=45)

#ifdef USE_MPI 
  ENDIF 
  n=LEN(value)
  CALL MPI_BCast(n,    1,MPI_INTEGER4,0,MPI_COMM_WORLD,ims_err) 
  CALL MPI_BCast(value,n,MPI_CHAR,0,MPI_COMM_WORLD,ims_err)
#endif 

  RETURN
END SUBROUTINE IO_READ_ASCII

! #######################################################################
! Write ASCII data; complete fields
! #######################################################################
SUBROUTINE IO_WRITE_ASCII_FIELD(fname, imax,jmax,kmax, u)
  
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro, ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  CHARACTER*(*),                    INTENT(IN) :: fname
  TINTEGER,                         INTENT(IN) :: imax,jmax,kmax
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(IN) :: u

! -----------------------------------------------------------------------
  TINTEGER idsp,kdsp, i,j,k
  CHARACTER*32 name_loc

! #######################################################################
#ifdef USE_MPI
  WRITE(name_loc,*) ims_pro; name_loc=TRIM(ADJUSTL(fname))//'-'//TRIM(ADJUSTL(name_loc))
  idsp = ims_offset_i; kdsp = ims_offset_k
#else
  name_loc = TRIM(ADJUSTL(fname))
  idsp = 0; kdsp = 0
#endif

  OPEN(unit=31,file=name_loc,status='unknown')

  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           WRITE(31,*) i+idsp,j,k+kdsp, u(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  CLOSE(31)

  RETURN
END SUBROUTINE IO_WRITE_ASCII_FIELD
