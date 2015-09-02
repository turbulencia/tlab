#include "types.h"

SUBROUTINE RAND_PDF(imax, jmax, kmax, seed, fsym, ipdf, u)

  IMPLICIT NONE

  TINTEGER imax, jmax, kmax
  TINTEGER seed, fsym, ipdf
  TREAL u(imax,jmax,kmax)

  TREAL RAN0, RANG
  TINTEGER ij,i,j,k,jmax_loc
  TREAL r0, r1

! defining symmetry or not
  IF ( fsym .EQ. 1 ) THEN; jmax_loc = jmax/2
  ELSE;                    jmax_loc = jmax;  ENDIF

! uniform pdf
  IF      ( ipdf .EQ. 1 ) THEN
     DO k=1,kmax
        DO ij=1,imax*jmax_loc
           u(ij,1,k) = RAN0(seed)-0.5
        ENDDO
     ENDDO

! gaussian pdf N(0,1)
  ELSE IF ( ipdf .EQ. 2 ) THEN
     r0 = C_0_R
     r1 = C_1_R
     DO k=1,kmax
        DO ij=1,imax*jmax_loc
           u(ij,1,k) = RANG(r0,r1,seed)
        ENDDO
     ENDDO
  ENDIF

! completing symmetric part if required
  IF ( fsym .EQ. 1 ) THEN
     DO k=1,kmax
        DO j=1,jmax/2
           DO i=1,imax
              u(i,jmax-j+1,k) = u(i,j,k)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE RAND_PDF
