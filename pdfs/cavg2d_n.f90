#include "types.h"
#include "dns_error.h"

#define NVARS_LOC 20
#define NMOMS_LOC 4
! the next macro must be the product of the previous two, for the format
#define FMT_LOC 80 

!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 2007/12/10 - J.P. Mellado
!#              Created
!# 2008/04/02 - J.P. Mellado
!#              Reformulation of the gate signal
!#
!########################################################################
!# DESCRIPTION
!#
!# Array b contains the conditioning field.
!# PDF format of TkStat file is used to allow variable b range at each plane.
!# The field gate gives a 2nd global conditioning (intermittency).
!# A last j-plane is added containing the CAVG constructed using all the
!# volume.
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

SUBROUTINE CAVG2D_N(fname, ibc, varname, igate, rtime, imax,jmax,kmax, &
     nvar, nbins, nmom, bmin, bmax, y, gate, b, data, cavg, wrk1d)

  USE DNS_TYPES,  ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname
  TREAL rtime
  TINTEGER imax,jmax,kmax, nvar, nbins, nmom, ibc
  TREAL y(jmax)
  TREAL b(imax,jmax,kmax), bmin, bmax
  TREAL cavg(nbins,nmom,nvar)
  TREAL wrk1d(nbins,*)
  CHARACTER*32 varname(nvar)

  INTEGER(1) gate(*), igate

  TYPE(pointers_dt), DIMENSION(nvar) :: data

! -------------------------------------------------------------------
  TINTEGER j, ip, iv, im
  TREAL xmin(NMOMS_LOC,NVARS_LOC), xmax(NMOMS_LOC,NVARS_LOC)

  CHARACTER*512 line1
  CHARACTER*32 str

#ifdef USE_MPI
  INTEGER ims_pro, ims_err

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

! ###################################################################
  IF ( NVARS_LOC .LT. nvar .OR. NMOMS_LOC .LT. nmom ) THEN
     CALL IO_WRITE_ASCII(efile, 'CAVG2D_N. Auxiliar array size too small')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! -------------------------------------------------------------------
! TkStat file; header
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
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
! Conditional average calculation of 1 variable in along planes
! ###################################################################
  DO j = 1,jmax

     DO iv = 1,nvar
        DO im = 1,nmom
           IF ( igate .GT. 0 ) THEN
              CALL CAVG2V2D1G(ibc, imax, jmax, kmax, j, igate, bmin, bmax, gate, b, &
                   data(iv)%field, nbins, im, wrk1d(1,1), cavg(1,im,iv), wrk1d(1,2),wrk1d(1,4))
           ELSE
              CALL CAVG2V2D(ibc, imax, jmax, kmax, j, bmin, bmax, b, &
                   data(iv)%field, nbins, im, wrk1d(1,1), cavg(1,im,iv), wrk1d(1,2),wrk1d(1,4))
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
        WRITE(21,1020) 1,j,ip,(cavg(ip,im,1),im=1,nmom*nvar),((xmin(im,iv),im=1,nmom),iv=1,nvar)
        DO ip=2, nbins-1
           WRITE(21,1010) 1,j, ip, (cavg(ip,im,1),im=1,nmom*nvar)
        ENDDO
        ip = nbins
        WRITE(21,1020) 1,j,ip,(cavg(ip,im,1),im=1,nmom*nvar),((xmax(im,iv),im=1,nmom),iv=1,nvar)

#ifdef USE_MPI
     ENDIF
#endif

  ENDDO

! ###################################################################
! Conditional average calculation of 1 variable in volume
! ###################################################################
  DO iv = 1,nvar
     DO im = 1,nmom
        IF ( igate .GT. 0 ) THEN
           CALL CAVG2V3D1G(ibc, imax, jmax, kmax, igate, bmin, bmax, gate, b, data(iv)%field, &
                nbins, im, wrk1d(1,1), cavg(1,im,iv), wrk1d(1,2), wrk1d(1,4))
        ELSE
           CALL CAVG2V3D(ibc, imax, jmax, kmax, bmin, bmax, b, data(iv)%field, &
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
     WRITE(21,1020) 1,jmax+1,ip,(cavg(ip,iv,1),iv=1,nmom*nvar),((xmin(im,iv),im=1,nmom),iv=1,nvar)
     DO ip=2, nbins-1
        WRITE(21,1010) 1, jmax+1, ip, (cavg(ip,iv,1),iv=1,nmom*nvar)
     ENDDO
     ip = nbins
     WRITE(21,1020) 1,jmax+1,ip,(cavg(ip,iv,1),iv=1,nmom*nvar),((xmax(im,iv),im=1,nmom),iv=1,nvar)

     CLOSE(21)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1010 FORMAT(I5,2(1X,I5),FMT_LOC(1X,E12.3E3))
1020 FORMAT(I5,2(1X,I5),FMT_LOC(1X,E12.3E3),FMT_LOC(1X,E12.3E3))

END SUBROUTINE CAVG2D_N


