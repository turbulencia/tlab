#include "types.h"
#include "dns_error.h"

!########################################################################
!# Chops string into list of strings
!########################################################################
SUBROUTINE LIST_STRING(line, n, a)
  IMPLICIT NONE

  CHARACTER*(*), INTENT(IN)    :: line
  TINTEGER,      INTENT(INOUT) :: n
  CHARACTER*(*), INTENT(OUT)   :: a(n)

  ! -------------------------------------------------------------------
  TINTEGER i, l1,l2,lmax

  ! ###################################################################
  i = 0                                             ! number of items

  lmax = LEN_TRIM(line)
  IF ( lmax > 0 ) THEN
    l1 = lmax -LEN_TRIM(ADJUSTL(line(1:lmax))) +1   ! absolute position of first nonblank in remaining string
    DO
      l2 = INDEX(line(l1:lmax),' ')                 ! relative position of first blank in remaining string

      i = i +1
      IF ( l2 == 0 ) THEN                           ! we found the last element
        a(i) = line(l1:lmax)
        EXIT
      ELSE
        a(i) = line(l1:l1+l2-1)
        l1 = lmax -LEN_TRIM(ADJUSTL(line(l1+l2:lmax))) +1
      END IF
      IF ( i == n ) EXIT

    END DO
  END IF

  n = i                                             ! return the number of items

  RETURN
END SUBROUTINE LIST_STRING

!########################################################################
!# Chops string into list of integers
!########################################################################
SUBROUTINE LIST_INTEGER(line, n, a)
  USE TLAB_CORE
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
  END IF
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
          END IF
          IF ( lloc .EQ. l2 ) THEN
            iread = .TRUE.
            ilast = lloc
          END IF
          IF ( iread ) THEN
            i = i + 1
            ! check the array is big enough
            IF ( i .GT. n ) THEN
              CALL TLAB_STOP(DNS_ERROR_PARAMETER)
            END IF
            READ(line(lfirst:ilast),*) a(i)
            lfirst = lloc
            GOTO 111
          END IF
        END DO
      END IF
111   IF ( lfirst .EQ. l2 ) GOTO 222
    END DO

    ! assign the correct size
222 n = i

    ! -------------------------------------------------------------------
    ! Matlab notation (first:step:last)
    ! -------------------------------------------------------------------
  ELSE
    READ(line(l1:lloc-1),*) a(1)
    l1 = lloc+1
    lloc = INDEX(line(l1:l2),':')
    IF ( lloc .EQ. 0 ) THEN
      CALL TLAB_STOP(DNS_ERROR_PARAMETER)
    END IF
    lloc = l1 + lloc - 1
    READ(line(l1:lloc-1),*) incr
    l1 = lloc+1
    READ(line(l1:l2),*) itmax

    ! check the array is big enough
    IF ( (itmax-a(1))/incr+1 .GT. n ) THEN
      CALL TLAB_STOP(DNS_ERROR_PARAMETER)
    ELSE
      n = (itmax-a(1))/incr+1
    END IF

    DO i = 2,n
      a(i) = a(1) + (i-1)*incr
    END DO

  END IF

  RETURN
END SUBROUTINE LIST_INTEGER

!########################################################################
!# Chops string into list of real numbers
!########################################################################
SUBROUTINE LIST_REAL(line, n, a)
  USE TLAB_CORE
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
  END IF
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
          END IF
          IF ( lloc .EQ. l2 ) THEN
            iread = .TRUE.
            ilast = lloc
          END IF
          IF ( iread ) THEN
            i = i + 1
            ! check the array is big enough
            IF ( i .GT. n ) THEN
              CALL TLAB_STOP(DNS_ERROR_PARAMETER)
            END IF
            READ(line(lfirst:ilast),*) a(i)
            lfirst = lloc
            GOTO 111
          END IF
        END DO
      END IF
111   IF ( lfirst .EQ. l2 ) GOTO 222
    END DO

    ! assign the correct size
222 n = i

    ! -------------------------------------------------------------------
    ! Matlab notation (first:step:last)
    ! -------------------------------------------------------------------
  ELSE
    READ(line(l1:lloc-1),*) a(1)
    l1 = lloc+1
    lloc = INDEX(line(l1:l2),':')
    IF ( lloc .EQ. 0 ) THEN
      CALL TLAB_STOP(DNS_ERROR_PARAMETER)
    END IF
    lloc = l1 + lloc - 1
    READ(line(l1:lloc-1),*) aincr
    l1 = lloc+1
    READ(line(l1:l2),*) amax


    ! check the array is big enough
    IF ( INT((amax-a(1))/aincr) +1 .GT. n ) THEN
      CALL TLAB_STOP(DNS_ERROR_PARAMETER)
    ELSE
      n = INT((amax-a(1))/aincr) +1
    END IF

    DO i = 2,n
      a(i) = a(1) + (i-1)*aincr
    END DO

  END IF

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
      END IF
    END DO
  END DO

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
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE SORT_REAL
