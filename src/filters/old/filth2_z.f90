#include "types.h"
#include "dns_error.h"

SUBROUTINE FILTH2_Z(iunifz, k1bc, imax,jmax,kmax, nz0,nz1, cf2z, z1, zf1, wrk)

  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER iunifz, k1bc, imax,jmax,kmax, nz0,nz1
  TREAL cf2z(*)
  TREAL z1(imax, jmax, kmax)
  TREAL zf1(imax, jmax, kmax)
  TREAL wrk(imax, jmax, *)

  TINTEGER nij, nz
  TINTEGER i2 

#ifdef USE_MPI
  TINTEGER npl
#endif

  i2 = 2
  nij = imax*jmax
  nz = nz0 + nz1

  IF ( MOD(nz0,i2) .NE. 0 .AND. MOD(nz1,i2) .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, 'FILTH2_Z. NZ2 is not even')
     CALL DNS_STOP(DNS_ERROR_LESEVEN)
  ENDIF

#ifdef USE_MPI
  IF ( ims_npro .GT. 1 ) THEN

! Trying 1-1 PE communication
     npl = nz/2
     IF ( npl .LE. kmax ) THEN
        CALL DNS_MPI_COPYPLN1(nij, kmax, npl, z1, wrk(1,1,1), wrk(1,1,1+npl+kmax))
     ELSE IF ( npl .LE. 2*kmax ) THEN
        CALL DNS_MPI_COPYPLN2(nij, kmax, npl, z1, wrk(1,1,1), wrk(1,1,1+npl+kmax))
     ELSE
        CALL IO_WRITE_ASCII(efile, 'FILTH2_Z. Size kmax too small for PARALLEL mode.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
        
     IF ( iunifz .EQ. 0 ) THEN
        IF ( k1bc .EQ. 0 ) THEN
           IF ( nz .EQ. 4 ) THEN
              CALL FILTH2MPPD4_MPI(kmax, nij, nz, cf2z, z1, zf1, wrk)
           ELSE IF ( nz .EQ. 6 ) THEN
              CALL FILTH2MPPD6_MPI(kmax, nij, nz, cf2z, z1, zf1, wrk)
           ELSE
              CALL FILTH2MPPD_MPI(kmax, nij, nz, cf2z, z1, zf1, wrk)
           ENDIF
        ENDIF
     ENDIF

  ELSE
#endif

     IF ( iunifz .EQ. 0 ) THEN
        IF ( k1bc .EQ. 0 ) THEN
           IF ( nz0+nz1 .EQ. 4 ) THEN
              CALL FILTH2MPPD4(kmax, nij, cf2z, z1, zf1)
           ELSE IF ( nz0+nz1 .EQ. 6 ) THEN
              CALL FILTH2MPPD6(kmax, nij, cf2z, z1, zf1)
           ELSE
              CALL FILTH2MPPD(kmax, nij, nz0, nz1, cf2z, z1, zf1)
           ENDIF
        ELSE
           CALL FILTH2MP(kmax, nij, nz0, nz1, cf2z, z1, zf1)
        ENDIF
     ELSE
        IF ( k1bc .EQ. 0 ) THEN
           CALL FILTH2MPPDNU(kmax, nij, nz0, nz1, cf2z, z1, zf1)
        ELSE
           CALL FILTH2MPNU(kmax, nij, nz0, nz1, cf2z, z1, zf1)
        ENDIF
     ENDIF

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE FILTH2_Z

