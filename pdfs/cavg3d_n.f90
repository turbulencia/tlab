#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 2008/04/02 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Simplified version of CAVG3D_N to consider only a 3D case
!#
!########################################################################
!# ARGUMENTS 
!#
!# igate    In     Gate level. If 0, no intermittency considered
!# nvar     In     Number of variables
!# nmom     In     Number of moments
!# ibc      In     BCs for array b: 1 for local bmin/bmax, 0 for fixed bmin/bmax 
!#
!########################################################################
SUBROUTINE CAVG3D_N(prefix, ibc, varname, igate, itime, rtime, imax, jmax, kmax, &
     nvar, nbins, nmom, bmin, bmax, gate, b, a, cavg, wrk1d)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

#define NVARS_LOC 20
#define NMOMS_LOC 4
! the next macro must be the product of the previous two, for the format
#define FMT_LOC 80 

  TINTEGER itime
  TREAL rtime
  TINTEGER imax, jmax, kmax, nvar, nbins, nmom, ibc
  TREAL b(imax,jmax,kmax), bmin, bmax
  TREAL a(imax,jmax,kmax,*)
  TREAL cavg(nbins,nmom,nvar)
  TREAL wrk1d(nbins,*)
  CHARACTER*32 prefix
  CHARACTER*32 varname(nvar)

  INTEGER(1) gate(*), igate

! -------------------------------------------------------------------
  TINTEGER ip, iv, im
  TREAL xmin(NMOMS_LOC,NVARS_LOC), xmax(NMOMS_LOC,NVARS_LOC)

  CHARACTER*32 fname, str
  CHARACTER*512 line1

#ifdef USE_MPI
  INTEGER ims_pro, ims_err

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

! ###################################################################
  IF ( NVARS_LOC .LT. nvar .OR. NMOMS_LOC .LT. nmom ) THEN
     CALL IO_WRITE_ASCII(efile, 'CAVG3D_N. Auxiliar array size too small')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! -------------------------------------------------------------------
! TkStat file; header
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     WRITE(fname,*) itime; fname=TRIM(ADJUSTL(prefix))//TRIM(ADJUSTL(fname))

     OPEN(unit=21,file=fname)

! comment section
     WRITE(21,'(A)') '# TkStat PDF format'
     IF ( igate .EQ. 0 ) THEN
        WRITE(21,'(A)') '# No conditioning on 3rd variable'
     ELSE
        WRITE(21,'(A37,I1)') &
             '# Conditioning on 3rd variable, level ', igate
     ENDIF

     WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime
     WRITE(21, '(A7,I8)') 'IMAX = ', i1
     WRITE(21, '(A7,I8)') 'JMAX = ', jmax+1

     line1 = 'I J N'
     DO iv = 1,nvar
        DO im = 1,nmom
           WRITE(str,*) im; line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(varname(iv)))//'Mom'//TRIM(ADJUSTL(str))
        ENDDO
     ENDDO
     DO iv = 1,nvar
        DO im = 1,nmom
           WRITE(str,*) im; line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(varname(iv)))//'Mom'//TRIM(ADJUSTL(str))//'_X'
        ENDDO
     ENDDO
     WRITE(21,'(A)') TRIM(ADJUSTL(line1))

#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
! Conditional average calculation of 1 variable in volume
! ###################################################################
  DO iv = 1,nvar
     DO im = 1,nmom
        IF ( igate .GT. 0 ) THEN
           CALL CAVG2V3D1G(ibc, imax, jmax, kmax, igate, bmin, bmax, gate, b, a(1,1,1,iv), &
                nbins, im, wrk1d(1,1), cavg(1,im,iv), wrk1d(1,2), wrk1d(1,4))
        ELSE
           CALL CAVG2V3D(ibc, imax, jmax, kmax, bmin, bmax, b, a(1,1,1,iv), &
                nbins, im, wrk1d(1,1), cavg(1,im,iv), wrk1d(1,2), wrk1d(1,4))
        ENDIF
        xmin(im,iv) = wrk1d(1,    1)
        xmax(im,iv) = wrk1d(nbins,1)

     ENDDO
  ENDDO

! -------------------------------------------------------------------
! TkStat output
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     ip = 1
     WRITE(21,1020) 1, 1, ip, (cavg(ip,iv,1),iv=1,nmom*nvar), (xmin(iv,1),iv=1,nmom*nvar)
     DO ip=2, nbins-1
        WRITE(21,1010) 1, 1, ip, (cavg(ip,iv,1),iv=1,nmom*nvar)
     ENDDO
     ip = nbins
     WRITE(21,1020) 1, 1, ip, (cavg(ip,iv,1),iv=1,nmom*nvar), (xmax(iv,1),iv=1,nmom*nvar)

     CLOSE(21)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1010 FORMAT(I5,2(1X,I5),FMT_LOC(1X,E12.3E3))
1020 FORMAT(I5,2(1X,I5),FMT_LOC(1X,E12.3E3),FMT_LOC(1X,E12.3E3))

END SUBROUTINE CAVG3D_N

