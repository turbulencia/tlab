#include "types.h"
#include "dns_error.h"
  
!########################################################################
!# Tool/Library UTILS
!#
!########################################################################
!# HISTORY
!#
!# 2000/02/01 - J.P. Mellado
!#              Created
!# 2007/02/25 - J.P. Mellado
!#              Modify name from CHOPI to LIST_INTEGER
!#
!########################################################################
!# DESCRIPTION
!#
!# Chops the string line containing integers into the elements of the 
!# array a
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE LIST_INTEGER(line, n, a)
  
  IMPLICIT NONE
  
  CHARACTER*(*),          INTENT(IN)    :: line
  TINTEGER,               INTENT(INOUT) :: n
  TINTEGER, DIMENSION(n), INTENT(OUT)   :: a
  
! -------------------------------------------------------------------
  TINTEGER i, lloc
  LOGICAL iread
  TINTEGER incr, itmax
  TINTEGER lfirst, ilast
  TINTEGER l1, l2
  
! ###################################################################
  l2 = LEN_TRIM(line)
  IF ( l2 .EQ. 0 ) THEN ! empty string
     n = 0
     RETURN
  ELSE
     l1 = l2 - LEN_TRIM(ADJUSTL(line)) + 1
  ENDIF
  lloc = INDEX(line(l1:l2),':')

! -------------------------------------------------------------------
! List separated by commas
! -------------------------------------------------------------------
  IF ( lloc .EQ. 0 ) THEN   
     i = 0
     lfirst = l1-1
     DO WHILE ( .TRUE. )
        lfirst = lfirst + 1
        IF ( (line(lfirst:lfirst) .NE. ' ' .AND. line(lfirst:lfirst) .NE. ',') ) THEN
! beggining of an item
           DO lloc = lfirst,l2
              iread = .FALSE.
              IF ( line(lloc:lloc) .EQ. ' ' .OR. line(lloc:lloc) .EQ. ',' ) THEN
                 iread = .TRUE.
                 ilast = lloc-1
              ENDIF
              IF ( lloc .EQ. l2 ) THEN
                 iread = .TRUE.
                 ilast = lloc
              ENDIF
              IF ( iread ) THEN
                 i = i + 1
! check the array is big enough
                 IF ( i .GT. n ) THEN
                    CALL DNS_STOP(DNS_ERROR_PARAMETER)
                 ENDIF
                 READ(line(lfirst:ilast),*) a(i)
                 lfirst = lloc
                 GOTO 111
              ENDIF
           ENDDO
        ENDIF
111     IF ( lfirst .EQ. l2 ) GOTO 222
     ENDDO
         
! assign the correct size
222  n = i
     
! -------------------------------------------------------------------
! Matlab notation (first:step:last)
! -------------------------------------------------------------------
  ELSE
     READ(line(l1:lloc-1),*) a(1)
     l1 = lloc+1
     lloc = INDEX(line(l1:l2),':')
     IF ( lloc .EQ. 0 ) THEN
        CALL DNS_STOP(DNS_ERROR_PARAMETER)
     ENDIF
     lloc = l1 + lloc - 1
     READ(line(l1:lloc-1),*) incr
     l1 = lloc+1
     READ(line(l1:l2),*) itmax
     
! check the array is big enough
     IF ( (itmax-a(1))/incr+1 .GT. n ) THEN
        CALL DNS_STOP(DNS_ERROR_PARAMETER)
     ELSE
        n = (itmax-a(1))/incr+1
     ENDIF
     
     DO i = 2,n
        a(i) = a(1) + (i-1)*incr
     ENDDO
     
  ENDIF
  
  RETURN
END SUBROUTINE LIST_INTEGER

!########################################################################
!########################################################################

SUBROUTINE LIST_REAL(line, n, a)

  IMPLICIT NONE

  CHARACTER*(*),          INTENT(IN)    :: line
  TINTEGER,               INTENT(INOUT) :: n
  TREAL, DIMENSION(n),    INTENT(OUT)   :: a

! -------------------------------------------------------------------
  TINTEGER i, lloc
  LOGICAL iread
  TREAL aincr, amax
  TINTEGER lfirst, ilast
  TINTEGER l1, l2

! ###################################################################
  l2 = LEN_TRIM(line)
  IF ( l2 .EQ. 0 ) THEN ! empty string
     n = 0
     RETURN
  ELSE
     l1 = l2 - LEN_TRIM(ADJUSTL(line)) + 1
  ENDIF
  lloc = INDEX(line(l1:l2),':')

! -------------------------------------------------------------------
! List separated by commas
! -------------------------------------------------------------------
  IF ( lloc .EQ. 0 ) THEN
     i = 0
     lfirst = l1-1
     DO WHILE ( .TRUE. )
        lfirst = lfirst + 1
        IF ( (line(lfirst:lfirst) .NE. ' ' .AND. line(lfirst:lfirst) .NE. ',') ) THEN
! beggining of an item
           DO lloc = lfirst,l2
              iread = .FALSE.
              IF ( line(lloc:lloc) .EQ. ' ' .OR. line(lloc:lloc) .EQ. ',' ) THEN
                 iread = .TRUE.
                 ilast = lloc-1
              ENDIF
              IF ( lloc .EQ. l2 ) THEN
                 iread = .TRUE.
                 ilast = lloc
              ENDIF
              IF ( iread ) THEN
                 i = i + 1
! check the array is big enough
                 IF ( i .GT. n ) THEN
                    CALL DNS_STOP(DNS_ERROR_PARAMETER)
                 ENDIF
                 READ(line(lfirst:ilast),*) a(i)
                 lfirst = lloc
                 GOTO 111
              ENDIF
           ENDDO
        ENDIF
111     IF ( lfirst .EQ. l2 ) GOTO 222
     ENDDO

! assign the correct size
222  n = i

! -------------------------------------------------------------------
! Matlab notation (first:step:last)
! -------------------------------------------------------------------
  ELSE
     READ(line(l1:lloc-1),*) a(1)
     l1 = lloc+1
     lloc = INDEX(line(l1:l2),':')
     IF ( lloc .EQ. 0 ) THEN
        CALL DNS_STOP(DNS_ERROR_PARAMETER)
     ENDIF
     lloc = l1 + lloc - 1
     READ(line(l1:lloc-1),*) aincr
     l1 = lloc+1
     READ(line(l1:l2),*) amax


! check the array is big enough
     IF ( INT((amax-a(1))/aincr) +1 .GT. n ) THEN
        CALL DNS_STOP(DNS_ERROR_PARAMETER)
     ELSE
        n = INT((amax-a(1))/aincr) +1
     ENDIF

     DO i = 2,n
        a(i) = a(1) + (i-1)*aincr
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE LIST_REAL

!########################################################################
!########################################################################

SUBROUTINE SORT_INTEGER(n,a) ! Sorting elements in array from min to max
  
  IMPLICIT NONE
  
  TINTEGER, INTENT(IN)                  :: n
  TINTEGER, DIMENSION(n), INTENT(INOUT) :: a

! -------------------------------------------------------------------
  TINTEGER i, j
  TINTEGER dummy
  
! ###################################################################
  DO i = 1,n-1
     DO j = i+1,n
        IF (a(j) .LT. a(i)) THEN
           dummy = a(i)
           a(i) = a(j)
           a(j) = dummy
        ENDIF
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE SORT_INTEGER

!########################################################################
!########################################################################

SUBROUTINE SORT_REAL(n,a)
  
  IMPLICIT NONE
  
  TINTEGER,            INTENT(IN)    :: n
  TREAL, DIMENSION(n), INTENT(INOUT) :: a
  
! -------------------------------------------------------------------
  TINTEGER i, j
  TREAL dummy
  
! ###################################################################
  DO i = 1,n-1
     DO j = i+1,n
        IF (a(j) .LT. a(i)) THEN
           dummy = a(i)
           a(i) = a(j)
           a(j) = dummy
        ENDIF
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE SORT_REAL
