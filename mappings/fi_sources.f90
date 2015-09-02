#include "types.h"

SUBROUTINE FI_SOURCES(flag_type, is, imax,jmax,kmax, y, s, wrk3d)

  USE DNS_GLOBAL, ONLY : scaley, ycoor_i

  IMPLICIT NONE

! -----------------------------------------------------------------------
  TINTEGER i,j,k
  TREAL xi, ycenter, yrel, dummy

! #######################################################################
  TINTEGER,                           INTENT(IN)    :: flag_type, is, imax,jmax,kmax
  TREAL, DIMENSION(*),                INTENT(IN)    :: y
  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN)    :: s
  TREAL, DIMENSION(imax*jmax*kmax),   INTENT(INOUT) :: wrk3d
  
  IF      ( flag_type .EQ. 1 ) THEN ! Left black to accomodate FI_RADIATION
  ELSE IF ( flag_type .EQ. 2 ) THEN ! Simple relaxation term
     ycenter=y(1)+scaley*ycoor_i(3)+0.005 !postition of s3 plus constant
     DO i=1,imax
        DO k=1,kmax
           DO j=1,jmax
              yrel=y(j)-ycenter
              xi=yrel/0.2 !thicknes constant
              wrk3d(i+(j-1)*imax+(k-1)*imax*jmax) = (1+TANH(xi))/2 !strength constant
           ENDDO
        ENDDO
     ENDDO
     dummy = C_1_R/0.02
     wrk3d(:) = - s(:,is)*wrk3d(:)*dummy
  ENDIF
  
  RETURN
END SUBROUTINE FI_SOURCES
