#include "types.h"
#include "dns_error.h"

#define NVARS_LOC 12

!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/27 - J.P. Mellado
!#              Created
!# 2007/09/27 - J.P. Mellado
!#              OpenDx added
!# 2008/01/17 - J.P. Mellado
!#              Conditional/Unconditional, not both 
!# 2008/04/03 - J.P. Mellado
!#              Gate reformulation
!#
!########################################################################
!# DESCRIPTION
!#
!# Array gate contains the global intermittency conditioning field.
!# Only the conditioned case is save for OpenDx analysis
!# A last j-plane is added containing the PDF constructed using all the
!# volume.
!# Note that amin/amax need to have required space !
!#
!########################################################################
!# ARGUMENTS 
!#
!# igate     In    Gate level. If 0, no intermittency considered
!# nvar      In    Number of variables
!# ibc       In    BCs: 0 homogeneous interval
!#                      1 local interval
!#                      2 local interval, analysis and drop no point
!#                      3 local interval, analysis and drop left point
!#                      4 local interval, analysis and drop right point
!#                      5 local interval, analysis and drop both points
!#
!########################################################################
SUBROUTINE PDF2D_N(fname, varname, igate, rtime, imax, jmax, kmax, &
     nvar, ibc, amin, amax, y, gate, data, nbins, npdf_size, pdf, wrk1d)

  USE DNS_TYPES,  ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname
  TREAL rtime
  TINTEGER imax, jmax, kmax, nvar
  TINTEGER ibc(nvar)
  TREAL amin(nvar), amax(nvar)
  TINTEGER nbins, npdf_size
  TREAL y(jmax)
  TREAL pdf(nbins,nvar)
  TREAL wrk1d(nbins,*) 
  CHARACTER*32 varname(nvar)

  INTEGER(1) gate(*), igate

  TYPE(pointers_dt), DIMENSION(nvar) :: data

! -------------------------------------------------------------------
  TINTEGER j, ip, nplim, iv, ibc_loc
  TREAL xmin(NVARS_LOC), xmax(NVARS_LOC), plim

  CHARACTER*512 line1

#ifdef USE_MPI
  INTEGER ims_pro, ims_err

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

! ###################################################################
  IF ( npdf_size .LT. nbins*nvar ) THEN
     CALL IO_WRITE_ASCII(efile, 'PDF2D_N. Working array size too small')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF
  IF ( NVARS_LOC .LT. nvar ) THEN
     CALL IO_WRITE_ASCII(efile, 'PDF2D_N. Aux array size too small')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! threshold in the PDF analysis
! This has to be reviewed for the case of small sample sizes, because
! maybe this value is never achieved !
  plim = C_1EM4_R

! -------------------------------------------------------------------
! TkStat file; header
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     OPEN(unit=21,file=fname)

! comment section
     WRITE(21,'(A)') '# PDF (normalization s.t. integral is 1)'
     IF ( igate .EQ. 0 ) THEN
        WRITE(21,'(A)') &
             '# No intermittency conditioning on 3rd variable'
     ELSE
        WRITE(21,'(A51,I1)') &
             '# Intermittency conditioning on gate, level ', igate
     ENDIF
     IF ( ibc(1) .GT. 1 ) THEN
        WRITE(21,'(A)') '# No PDF analysis'
     ELSE
        WRITE(21,'(A26,E14.7E3)') '# PDF analysis, threshold ', plim
     ENDIF

     WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime
     WRITE(21, '(A7,I8)') 'IMAX = ', i1
     WRITE(21, '(A7,I8)') 'JMAX = ', jmax+1

     line1 = 'I J IP'
     DO iv = 1,nvar
        line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(varname(iv)))
     ENDDO
     DO iv = 1,nvar
        line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(varname(iv)))//'_X'
     ENDDO
     WRITE(21,'(A)') TRIM(ADJUSTL(line1))

#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
! PDF calculation of 1 variable along planes
! ###################################################################
  DO j = 1,jmax

     DO iv = 1,nvar
        IF ( igate .EQ. 0 ) THEN
           CALL PDF1V2D(i1, ibc(iv), imax, jmax, kmax, j, amin(iv), amax(iv), &
                data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
        ELSE
           CALL PDF1V2D1G(i1, ibc(iv), imax, jmax, kmax, j, igate, amin(iv), amax(iv), &
                gate, data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
        ENDIF

! threshold for analysis set s.t. single points are removed
        IF ( ibc(iv) .GT. 1 ) THEN
           ibc_loc = ibc(iv)-2
           CALL PDF_ANALIZE&
                (nbins, ibc_loc, xmin(iv), xmax(iv), pdf(1,iv), plim, amin(iv), amax(iv), nplim)
           IF ( igate .EQ. 0 ) THEN
              CALL PDF1V2D(i1, i0, imax, jmax, kmax, j, amin(iv), amax(iv), &
                   data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
           ELSE
              CALL PDF1V2D1G(i1, i0, imax, jmax, kmax, j, igate, amin(iv), amax(iv), &
                   gate, data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
           ENDIF
        ENDIF

     ENDDO

! -------------------------------------------------------------------
! TkStat output
! -------------------------------------------------------------------
#ifdef USE_MPI
     IF ( ims_pro .EQ. 0 ) THEN
#endif
        ip = 1
        WRITE(21,1020) 1, j, ip, (pdf(ip,iv),iv=1,nvar), (xmin(iv),iv=1,nvar)
        DO ip=2, nbins-1
           WRITE(21,1010) 1, j, ip, (pdf(ip,iv),iv=1,nvar)
        ENDDO
        ip = nbins
        WRITE(21,1020) 1, j, ip, (pdf(ip,iv),iv=1,nvar), (xmax(iv),iv=1,nvar)

#ifdef USE_MPI
     ENDIF
#endif
  ENDDO

! ###################################################################
! PDF calculation of 1 variable in 3D space
! ###################################################################
  DO iv = 1,nvar
     IF ( igate .EQ. 0 ) THEN
        CALL PDF1V3D(i1, ibc(iv), imax, jmax, kmax, amin(iv), amax(iv), &
             data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
     ELSE
        CALL PDF1V3D1G(i1, ibc(iv), imax, jmax, kmax, igate, amin(iv), amax(iv), &
             gate, data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
     ENDIF

! threshold for analysis set s.t. single points are removed
     IF ( ibc(iv) .GT. 1 ) THEN
        ibc_loc = ibc(iv)-2
        CALL PDF_ANALIZE&
             (nbins, ibc_loc, xmin(iv), xmax(iv), pdf(1,iv), plim, amin(iv), amax(iv), nplim)
        IF ( igate .EQ. 0 ) THEN
           CALL PDF1V3D(i1, i0, imax, jmax, kmax, amin(iv), amax(iv), &
                data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
        ELSE
           CALL PDF1V3D1G(i1, i0, imax, jmax, kmax, igate, amin(iv), amax(iv), &
                gate, data(iv)%field, nbins, xmin(iv), xmax(iv), pdf(1,iv), wrk1d)
        ENDIF
     ENDIF
  ENDDO

! -------------------------------------------------------------------
! TkStat output
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     ip = 1
     WRITE(21,1020) 1, jmax+1, ip, (pdf(ip,iv),iv=1,nvar), (xmin(iv),iv=1,nvar)
     DO ip=2, nbins-1
        WRITE(21,1010) 1, jmax+1, ip, (pdf(ip,iv),iv=1,nvar)
     ENDDO
     ip = nbins
     WRITE(21,1020) 1, jmax+1, ip, (pdf(ip,iv),iv=1,nvar), (xmax(iv),iv=1,nvar)

     CLOSE(21)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1010 FORMAT(I5,2(1X,I5),NVARS_LOC(1X,E12.3E3))
1020 FORMAT(I5,2(1X,I5),NVARS_LOC(1X,E12.3E3),NVARS_LOC(1X,E12.3E3))

END SUBROUTINE PDF2D_N


