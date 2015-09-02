#include "types.h"

!########################################################################
!# Tool/Library 
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/04 - J.P. Mellado
!#              Created
!# 2012/09/30 - J.P. Mellado
!#              Adding INT1 routines
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, amin, y, a, at, surface, wrk2d)

  IMPLICIT NONE

  TINTEGER imax, jmax, kmax, jmin_loc
  TREAL y(jmax), amin
  TREAL a(*), at(jmax,imax*kmax)
  TREAL surface(*), wrk2d(imax*kmax)

! -------------------------------------------------------------------
  TINTEGER i, j, jkmax, ikmax

! ###################################################################
  jkmax = jmax*kmax
  ikmax = imax*kmax

! -------------------------------------------------------------------
! Make y direction the first one; x becomes the last one
! -------------------------------------------------------------------
  CALL DNS_TRANSPOSE(a, imax, jkmax, imax, at, jkmax)

  DO i = 1,ikmax
     DO j = jmin_loc+1,jmax
        IF ( at(j,i) .GT. amin ) EXIT
     ENDDO
     IF ( j .EQ. jmax+1 ) THEN
        wrk2d(i) = y(jmax)
     ELSE
        wrk2d(i) = y(j-1) + (y(j)-y(j-1))/(at(j,i)-at(j-1,i))*(amin-at(j-1,i))
     ENDIF
  ENDDO

! -------------------------------------------------------------------
! Put array in right order
! -------------------------------------------------------------------
  CALL DNS_TRANSPOSE(wrk2d, kmax, imax, kmax, surface, imax)

  RETURN
END SUBROUTINE SL_LOWER_BOUNDARY

!########################################################################
!########################################################################

SUBROUTINE SL_UPPER_BOUNDARY(imax,jmax,kmax, jmax_loc, amin, y, a, at, surface, wrk2d)
  
  IMPLICIT NONE

  TINTEGER imax, jmax, kmax, jmax_loc
  TREAL y(jmax), amin
  TREAL a(*), at(jmax,imax*kmax)
  TREAL surface(*), wrk2d(imax*kmax)

! -------------------------------------------------------------------
  TINTEGER i, j, jkmax, ikmax

! ###################################################################
  jkmax = jmax*kmax
  ikmax = imax*kmax

! -------------------------------------------------------------------
! Make y direction the first one; x becomes the last one
! -------------------------------------------------------------------
  CALL DNS_TRANSPOSE(a, imax, jkmax, imax, at, jkmax)

  DO i = 1,ikmax
     DO j = jmax_loc-1,1,-1
        IF ( at(j,i) .GT. amin ) EXIT
     ENDDO
     IF ( j .EQ. 0 ) THEN
        wrk2d(i) = y(1)
     ELSE
        wrk2d(i) = y(j+1) + (y(j)-y(j+1))/(at(j,i)-at(j+1,i))*(amin-at(j+1,i))
     ENDIF
  ENDDO

! -------------------------------------------------------------------
! Put array in right order
! -------------------------------------------------------------------
  CALL DNS_TRANSPOSE(wrk2d, kmax, imax, kmax, surface, imax)

  RETURN
END SUBROUTINE SL_UPPER_BOUNDARY

!########################################################################
!########################################################################
SUBROUTINE BOUNDARY_LOWER_INT1(imax,jmax,kmax, avalue, y, a, at, surface, wrk2d)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax
  TREAL,      DIMENSION(jmax),           INTENT(IN)    :: y
  INTEGER(1),                            INTENT(IN)    :: avalue
  INTEGER(1), DIMENSION(jmax,imax*kmax), INTENT(IN)    :: a ! real shape is imax*jmax*kmax
  INTEGER(1), DIMENSION(jmax,imax*kmax), INTENT(INOUT) :: at(jmax,imax*kmax)
  TREAL,      DIMENSION(imax*kmax),      INTENT(INOUT) :: wrk2d
  TREAL,      DIMENSION(imax*kmax),      INTENT(OUT)   :: surface

! -------------------------------------------------------------------
  TINTEGER i, j, jkmax, ikmax

! ###################################################################
  jkmax = jmax*kmax
  ikmax = imax*kmax

! Make y direction the first one; x becomes the last one
  IF ( imax .GT. 1 ) THEN
     CALL DNS_TRANSPOSE_INT1(a, imax, jkmax, imax, at, jkmax)
  ELSE
     at = a
  ENDIF

  DO i = 1,ikmax
     DO j = 1,jmax
        IF ( at(j,i) .EQ. avalue ) EXIT
     ENDDO
     IF ( j .EQ. jmax+1 ) THEN
        wrk2d(i) = y(jmax) ! in case the avalue is never attained
     ELSE
        wrk2d(i) = y(j)
     ENDIF
  ENDDO

! Put array in right order
  CALL DNS_TRANSPOSE(wrk2d, kmax, imax, kmax, surface, imax)

  RETURN
END SUBROUTINE BOUNDARY_LOWER_INT1

! ###################################################################
! ###################################################################

SUBROUTINE BOUNDARY_UPPER_INT1(imax,jmax,kmax, avalue, y, a, at, surface, wrk2d)
  
  IMPLICIT NONE

  TINTEGER imax,jmax,kmax
  TREAL,      DIMENSION(jmax),           INTENT(IN)    :: y
  INTEGER(1),                            INTENT(IN)    :: avalue
  INTEGER(1), DIMENSION(jmax,imax*kmax), INTENT(IN)    :: a ! real shape is imax*jmax*kmax
  INTEGER(1), DIMENSION(jmax,imax*kmax), INTENT(INOUT) :: at(jmax,imax*kmax)
  TREAL,      DIMENSION(imax*kmax),      INTENT(INOUT) :: wrk2d
  TREAL,      DIMENSION(imax*kmax),      INTENT(OUT)   :: surface

! -------------------------------------------------------------------
  TINTEGER i, j, jkmax, ikmax

! ###################################################################
  jkmax = jmax*kmax
  ikmax = imax*kmax

! Make y direction the first one; x becomes the last one
  IF ( imax .GT. 1 ) THEN
     CALL DNS_TRANSPOSE_INT1(a, imax, jkmax, imax, at, jkmax)
  ELSE
     at = a
  ENDIF

  DO i = 1,ikmax
     DO j = jmax,1,-1
        IF ( at(j,i) .EQ. avalue ) EXIT
     ENDDO
     IF ( j .EQ. 0 ) THEN; wrk2d(i) = y(  1) ! in case the avalue is never attained
     ELSE;                 wrk2d(i) = y(MIN(j+1,jmax)); ENDIF ! +1 so that it coincides with lower boundary
  ENDDO

! Put array in right order
  CALL DNS_TRANSPOSE(wrk2d, kmax, imax, kmax, surface, imax)

  RETURN
END SUBROUTINE BOUNDARY_UPPER_INT1

! ###################################################################
! ###################################################################

FUNCTION UPPER_THRESHOLD(jmax, uc, u, y)
  
  IMPLICIT NONE
  
  TINTEGER,               INTENT(IN) :: jmax
  TREAL, DIMENSION(jmax), INTENT(IN) :: u, y
  TREAL,                  INTENT(IN) :: uc
  TREAL UPPER_THRESHOLD
  
! -------------------------------------------------------------------
  TINTEGER j
  
! ###################################################################
  UPPER_THRESHOLD = C_0_R
  DO j = jmax-1,1,-1
     IF ( (u(j)-uc)*(u(j+1)-uc) .LT. C_0_R ) THEN
        UPPER_THRESHOLD = y(j) + (uc-u(j))*(y(j+1)-y(j))/(u(j+1)-u(j))
        EXIT
     ENDIF
  ENDDO
  
  RETURN
END FUNCTION UPPER_THRESHOLD

! ###################################################################
! ###################################################################

FUNCTION LOWER_THRESHOLD(jmax, uc, u, y)
  
  IMPLICIT NONE
  
  TINTEGER,               INTENT(IN) :: jmax
  TREAL, DIMENSION(jmax), INTENT(IN) :: u, y
  TREAL,                  INTENT(IN) :: uc
  TREAL LOWER_THRESHOLD
  
! -------------------------------------------------------------------
  TINTEGER j
  
! ###################################################################
  LOWER_THRESHOLD = C_0_R
  DO j = 1,jmax-1
     IF ( (u(j)-uc)*(u(j+1)-uc) .LT. C_0_R ) THEN
        LOWER_THRESHOLD = y(j) + (uc-u(j))*(y(j+1)-y(j))/(u(j+1)-u(j))
        EXIT
     ENDIF
  ENDDO
  
  RETURN
END FUNCTION LOWER_THRESHOLD

! ###################################################################
! ###################################################################

SUBROUTINE DELTA_X(imax, jmax, y, a, delta, delta_d, delta_u, A2, eta)

  IMPLICIT NONE

  TINTEGER,                    INTENT(IN)  :: imax, jmax
  TREAL, DIMENSION(jmax),      INTENT(IN)  :: y
  TREAL, DIMENSION(imax,jmax), INTENT(IN)  :: a
  TREAL,                       INTENT(IN)  :: A2, eta
  TREAL, DIMENSION(imax),      INTENT(OUT) :: delta_d, delta_u, delta

! -------------------------------------------------------------------
  TINTEGER i, j
  TREAL DA, A_05, y_center

! ###################################################################
  y_center = C_05_R*(y(jmax/2)+y(jmax/2+1))

  DO i = 1,imax
     DA = C_05_R*(a(i,jmax/2)+ a(i,jmax/2+1)) - A2
     A_05 = A2 + eta*DA

     DO j = 1,jmax/2
        IF ( a(i,j  ) .LE. A_05 .AND. a(i,j+1) .GT. A_05 ) THEN
           delta_d(i) = y(j) + (A_05-a(i,j))*(y(j+1)-y(j))/(a(i,j+1)-a(i,j))
        ENDIF
     ENDDO
     delta_d(i) = y_center - delta_d(i)

     DO j = jmax/2+1,jmax
        IF ( a(i,j  ) .GT. A_05 .AND. a(i,j+1) .LE. A_05 ) THEN
           delta_u(i) = y(j) + (A_05-a(i,j))*(y(j+1)-y(j))/(a(i,j+1)-a(i,j))
        ENDIF
     ENDDO
     delta_u(i) = delta_u(i) - y_center

     delta(i) = C_05_R*(delta_d(i)+delta_u(i))
  ENDDO

  RETURN
END SUBROUTINE DELTA_X
