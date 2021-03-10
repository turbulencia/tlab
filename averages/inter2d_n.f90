#include "types.h"
#include "dns_error.h"

#define NVARS_LOC 20

!########################################################################
!#
!# Intermittency factors, i.e., area fraction occupied by each gate level
!#
!########################################################################
SUBROUTINE INTER2D_N(fname, parname, rtime, nx,ny,nz, npar, y, gate)

USE DNS_CONSTANTS, ONLY : efile, lfile

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname, parname(npar)
  TREAL rtime
  TINTEGER,           INTENT(IN)    :: nx,ny,nz, npar ! npar is the number of partitions in gate field
  TREAL,              INTENT(IN)    :: y(ny)          ! heights of each plane
  INTEGER(1),         INTENT(IN)    :: gate(*)        ! field with partitions

  ! -------------------------------------------------------------------
  TINTEGER ip, j
  TREAL inter(NVARS_LOC), INTER1V2D
  INTEGER(1) gate_level

  CHARACTER*512 line1

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  IF ( NVARS_LOC .LT. npar ) THEN
    CALL IO_WRITE_ASCII(efile, 'INTER2D_N. Aux array size too small')
    CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

  ! -------------------------------------------------------------------
  ! TkStat file; header
  ! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
    OPEN(unit=21,file=fname)

    WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime
    WRITE(21, '(A7,I8)') 'IMAX = ', 1
    WRITE(21, '(A7,I8)') 'JMAX = ', ny

    line1 = 'I J Y'
    DO ip = 1,npar
      line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(parname(ip)))
    ENDDO
    WRITE(21,'(A)') TRIM(ADJUSTL(line1))

#ifdef USE_MPI
  ENDIF
#endif

  ! ###################################################################
  DO j = 1,ny

    DO ip = 1,npar
      gate_level = INT(ip,KIND=1)
      inter(ip) = INTER1V2D(nx,ny,nz, j, gate_level, gate)
    ENDDO

#ifdef USE_MPI
    IF ( ims_pro .EQ. 0 ) THEN
#endif
      WRITE(21,1010) 1, j, y(j), (inter(ip),ip=1,npar)
#ifdef USE_MPI
    ENDIF
#endif

  ENDDO

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
    CLOSE(21)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1010 FORMAT(I5,(1X,I5),(1X,G_FORMAT_R),NVARS_LOC(1X,G_FORMAT_R))

END SUBROUTINE INTER2D_N
