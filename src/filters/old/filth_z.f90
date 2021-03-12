#include "types.h"
#include "dns_error.h"

SUBROUTINE FILTH_Z(iunifz, k1bc, imax,jmax,kmax, nz, cfz, z1, zf1, wrk)
  
  USE DNS_CONSTANTS, ONLY : efile

#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif


  TINTEGER iunifz, k1bc, imax,jmax,kmax, nz
  TREAL cfz(*)
  TREAL z1(imax, jmax, kmax)
  TREAL zf1(imax, jmax, kmax)
  TREAL wrk(imax, jmax, *)

  TINTEGER nij
  TINTEGER i2 

#ifdef USE_MPI
  TINTEGER npl
#endif

  i2 = 2
  nij = imax*jmax

  IF ( MOD(nz,i2) .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, 'FILTH_Z. NZ is not even')
     CALL DNS_STOP(DNS_ERROR_LESEVEN)
  ENDIF

#ifdef USE_MPI
  IF ( ims_npro .GT. 1 ) THEN

! 1-1 PE communication
     npl = nz/2
     IF ( npl .LE. kmax ) THEN
        CALL DNS_MPI_COPYPLN1(nij, kmax, npl, z1, wrk(1,1,1), wrk(1,1,1+npl+kmax))
     ELSE IF ( npl .LE. 2*kmax ) THEN
        CALL DNS_MPI_COPYPLN2(nij, kmax, npl, z1, wrk(1,1,1), wrk(1,1,1+npl+kmax))
     ELSE
        CALL IO_WRITE_ASCII(efile, 'FILTH_Z. Size kmax too small for PARALLEL mode.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

     IF ( iunifz .EQ. 0 ) THEN
        IF ( k1bc .EQ. 0 ) THEN
           IF ( nz .EQ. 2 ) THEN
              CALL FILTHMPPD2_MPI(kmax, nij, nz, z1, zf1, wrk)
           ELSE IF ( nz .EQ. 4 ) THEN
              CALL FILTHMPPD4_MPI(kmax, nij, nz, z1, zf1, wrk)
           ELSE
              CALL FILTHMPPD_MPI(kmax, nij, nz, z1, zf1, wrk)
           ENDIF
        ENDIF
     ENDIF

  ELSE
#endif

     IF ( iunifz .EQ. 0 ) THEN
        IF ( k1bc .EQ. 0 ) THEN
           IF ( nz .EQ. 2 ) THEN
              CALL FILTHMPPD2(kmax, nij, z1, zf1)
           ELSE IF ( nz .EQ. 4 ) THEN
              CALL FILTHMPPD4(kmax, nij, z1, zf1)
           ELSE
              CALL FILTHMPPD(kmax, nij, nz, z1, zf1)
           ENDIF
        ELSE
           CALL FILTHMP(kmax, nij, nz, cfz, z1, zf1)
        ENDIF
     ELSE
        IF ( k1bc .EQ. 0 ) THEN
           CALL FILTHMPPDNU(kmax, nij, nz, cfz, z1, zf1)
        ELSE
           IF ( nz .EQ. 2 ) THEN
              CALL FILTHMPNU2(kmax, nij, cfz, z1, zf1)
           ELSE IF ( nz .EQ. 4 ) THEN
              CALL FILTHMPNU4(kmax, nij, cfz, z1, zf1)
           ELSE IF ( nz .EQ. 6 ) THEN
              CALL FILTHMPNU6(kmax, nij, cfz, z1, zf1)
           ELSE
              CALL FILTHMPNU(kmax, nij, nz, cfz, z1, zf1)
           ENDIF
        ENDIF
     ENDIF

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE FILTH_Z

