SUBROUTINE REYFLUCT2D(imax, jmax, kmax, dx, dz, area, u)

  IMPLICIT NONE

#include "types.h"

  TINTEGER imax, jmax, kmax
  TREAL dx(imax)
  TREAL dz(kmax)
  TREAL area
  TREAL u(imax,jmax,kmax)
  TREAL umn
  TINTEGER i,j,k
  TREAL AVG_IK

  DO j=1, jmax
     umn = AVG_IK(imax, jmax, kmax, j, u, dx, dz, area)
     DO k=1, kmax
        DO i=1,imax
           u(i,j,k) = u(i,j,k) - umn
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REYFLUCT2D
