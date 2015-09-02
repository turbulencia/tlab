#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/12/04 - J.P. Mellado
!#              Created
!# 2008/04/02 - J.P. Mellado
!#              Reformulation of the gate signal
!#
!########################################################################
!# DESCRIPTION
!#
!# Derived from CPDFFILE to consider in general N variables passed 
!# through array a(imax,jmax,kmax,nvar). Array b contains the conditioning 
!# field.
!# 
!########################################################################
!# ARGUMENTS 
!#
!# igate    In     Gate level. If 0, no intermittency considered
!# nvar     In     Number of variables
!#
!########################################################################
SUBROUTINE CPDF3D_N(fname, varname, igate, inorm, ianalyze, rtime, &
     imax, jmax, kmax, nvar, npz, nbins, gate, b, a, pdf, wrk1d)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif
  TINTEGER inorm, ianalyze
  TINTEGER imax, jmax, kmax, nvar
  TINTEGER nbins, npz
  TREAL rtime
  TREAL b(*)
  TREAL a(imax,jmax,kmax,*)
  TREAL pdf(nbins,nvar)
  TREAL wrk1d(nbins,3)

  CHARACTER*32 fname
  CHARACTER*32 varname(nvar)

  INTEGER(1) gate(*), igate

! -------------------------------------------------------------------
  TREAL xmin(nvar), xmax(nvar)
  TINTEGER ip, nplim, iv, nsample
  TREAL amin, amax, plim

  TINTEGER iz
  TREAL bmin, bmax, bdelta, bcenter

  CHARACTER*512 line1

#ifdef USE_MPI
  INTEGER ims_pro, ims_err

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

#define L_FORMAT_MAX 20

! ###################################################################
  IF ( nvar .GT. L_FORMAT_MAX ) THEN
     CALL IO_WRITE_ASCII(efile, 'CPDF3D_N. Format length too short.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! threshold in the PDF analysis
! This has to be reviewed for the case of small sample sizes, because
! maybe this value is never achieved !
  plim = C_1EM3_R

! calculate min/max of conditioning scalar
  CALL MINMAX(imax,jmax,kmax, b, bmin,bmax)
  bdelta = (bmax-bmin)/M_REAL(npz)

! eliminate the two outer bins
!  bmin = bmin + bdelta
!  bmax = bmax - bdelta
!  bdelta = (bmax-bmin)/M_REAL(npz)

! -------------------------------------------------------------------
! Initialize TkStat output
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     OPEN(unit=21,file=fname)

! comment section
     WRITE(21,'(A23,2(1X,E14.7E3))') '# Conditioning between ',bmin,bmax

! not yet implemented
     IF ( inorm .EQ. 0 ) THEN
        WRITE(21,'(A)') '# Histogram (no normalization)'
     ELSE
        WRITE(21,'(A)') '# PDF (normalization s.t. integral is 1)'
     ENDIF
     IF ( igate .EQ. 0 ) THEN
        WRITE(21,'(A)') &
             '# No intermittency conditioning on 3rd variable'
     ELSE
        WRITE(21,'(A51,I1)') &
             '# Intermittency conditioning on gate, level ', igate
     ENDIF
     IF ( ianalyze .EQ. 0 ) THEN
        WRITE(21,'(A)') '# No PDF analysis'
     ELSE
        WRITE(21,'(A26,E14.7E3)') '# PDF analysis, threshold ', plim
     ENDIF

     WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime

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
! PDF calculation of 1 variable in 3D space
! ###################################################################
  DO iz = 1, npz
     bcenter = bmin + bdelta*C_05_R + (iz-1)*bdelta

! loop on all variables
     DO iv = 1, nvar
        IF ( igate .GT. 0 ) THEN
           CALL CPDF2V3D1G(inorm, i1, imax, jmax, kmax, igate, amin, amax, gate, b, a(1,1,1,iv),&
                bcenter, bdelta, nbins, wrk1d(1,1), pdf(1,iv), wrk1d(1,2), nsample)
        ELSE
           CALL CPDF2V3D(inorm, i1, imax, jmax, kmax, amin, amax, b, a(1,1,1,iv),&
                bcenter, bdelta, nbins, wrk1d(1,1), pdf(1,iv), wrk1d(1,2), nsample)
        ENDIF
        xmin(iv) = wrk1d(1,    1)
        xmax(iv) = wrk1d(nbins,1)

! threshold for analysis set s.t. single points are removed
        IF ( ianalyze .EQ. 1 ) THEN
           CALL PDF_ANALIZE(nbins, i0, xmin(iv), xmax(iv), pdf(1,iv), plim, amin, amax, nplim)
           IF ( igate .GT. 0 ) THEN
              CALL CPDF2V3D1G(inorm, i0, imax,jmax,kmax, igate, amin, amax, gate, b, a(1,1,1,iv),&
                   bcenter, bdelta, nbins, wrk1d(1,1), pdf(1,iv), wrk1d(1,2), nsample)
           ELSE
              CALL CPDF2V3D(inorm, i0, imax, jmax, kmax, amin, amax, b, a(1,1,1,iv),&
                   bcenter, bdelta, nbins, wrk1d(1,1), pdf(1,iv), wrk1d(1,2), nsample)
           ENDIF
           xmin(iv) = wrk1d(1,    1)
           xmax(iv) = wrk1d(nbins,1)
        ENDIF

     ENDDO

! -------------------------------------------------------------------
! TkStat output; body
! -------------------------------------------------------------------
#ifdef USE_MPI
     IF ( ims_pro .EQ. 0 ) THEN
#endif
        ip = 1
        WRITE(21,1030) 1, iz, ip, (pdf(ip,iv), iv=1,nvar), (xmin(iv),iv=1,nvar)

        DO ip=2, nbins-1
           WRITE(21,1020) 1, iz, ip, (pdf(ip,iv), iv=1,nvar)
        ENDDO

        ip = nbins
        WRITE(21,1030) 1, iz, ip, (pdf(ip,iv), iv=1,nvar), (xmax(iv),iv=1,nvar)

1020    FORMAT(I5,2(1X,I5),L_FORMAT_MAX(1X,E12.3E3))
1030    FORMAT(I5,2(1X,I5),L_FORMAT_MAX(1X,E12.3E3),L_FORMAT_MAX(1X,E12.3E3))

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
END SUBROUTINE CPDF3D_N
