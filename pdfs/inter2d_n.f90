#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 2008/01/17 - J.P. Mellado
!#              Created
!# 2008/04/02 - J.P. Mellado
!#              Reformulation of the gate signal
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# npar    In    Number of partitions of the gate field
!#
!########################################################################
SUBROUTINE INTER2D_N(fname, parname, rtime, imax,jmax,kmax, npar, y, gate, gate_levels)

#ifdef USE_MPI
  USE DNS_MPI
#endif
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

#define NVARS_LOC 20

  CHARACTER*(*) fname
  TREAL rtime
  TINTEGER imax, jmax, kmax, npar
  TREAL y(jmax)
  CHARACTER*32 parname(npar)

  INTEGER(1) gate(*), gate_levels(npar)

! -------------------------------------------------------------------
  TINTEGER ip, j
  TREAL inter(NVARS_LOC)

  CHARACTER*512 line1

! ###################################################################
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
     WRITE(21, '(A7,I8)') 'IMAX = ', i1
     WRITE(21, '(A7,I8)') 'JMAX = ', jmax

     line1 = 'I J Y'
     DO ip = 1,npar
        line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(parname(ip)))
     ENDDO
     WRITE(21,'(A)') TRIM(ADJUSTL(line1))

#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
  DO j = 1,jmax

     DO ip = 1,npar
        CALL INTER1V2D(imax, jmax, kmax, j, gate_levels(ip), gate, inter(ip))
     ENDDO

! -------------------------------------------------------------------
! TkStat output
! -------------------------------------------------------------------
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

1010 FORMAT(I5,(1X,I5),(1X,E12.5E3),NVARS_LOC(1X,E12.5E3))

  RETURN
END SUBROUTINE INTER2D_N
