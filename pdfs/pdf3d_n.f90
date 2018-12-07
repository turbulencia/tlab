#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/19 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Simplified version of PDF2D_N to consider only a 3D case
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE PDF3D_N(fname, varname, ianalyze, rtime, &
     imax, jmax, kmax, nvar, nbins, a, pdf, wrk1d)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ianalyze
  TREAL rtime
  TINTEGER imax, jmax, kmax, nvar
  TINTEGER nbins
  TREAL a(imax,jmax,kmax,nvar)
  TREAL pdf(nbins+2,nvar)
  TREAL wrk1d(*)
  CHARACTER*32 fname
  CHARACTER*32 varname(nvar)

! -------------------------------------------------------------------
  TINTEGER ip, nplim, iv
  TREAL amin, amax, plim

  CHARACTER*512 line1
#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

#define L_FORMAT_MAX 20

! ###################################################################
  IF ( nvar .GT. L_FORMAT_MAX ) THEN
     CALL IO_WRITE_ASCII(efile, 'PDF3D_N. Format length too short.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! threshold in the PDF analysis
! This has to be reviewed for the case of small sample sizes, because
! maybe this value is never achieved !
  plim = C_1EM4_R

! ###################################################################
! PDF calculation of 1 variable in 3D space
! ###################################################################
  DO iv = 1, nvar
     CALL PDF1V3D(i1, imax, jmax, kmax, &
          amin, amax, a(1,1,1,iv), nbins, pdf(1,iv), wrk1d)

! threshold for analysis set s.t. single points are removed
     IF ( ianalyze .EQ. 1 ) THEN
        CALL PDF_ANALIZE(nbins, i0, pdf(1,iv), plim, amin, amax, nplim)
        CALL PDF1V3D(i0, imax, jmax, kmax, &
             amin, amax, a(1,1,1,iv), nbins, pdf(1,iv), wrk1d)
     ENDIF
  ENDDO

! ###################################################################
! TkStat output
! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

     OPEN(unit=21,file=fname)

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
! comment section
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

! -------------------------------------------------------------------
! body
! -------------------------------------------------------------------
     ip = 1
     WRITE(21,1030) 1, 1, ip, (pdf(ip,iv), iv=1,nvar), (pdf(nbins+1,iv),iv=1,nvar)

     DO ip = 2,nbins-1
        WRITE(21,1020) 1, 1, ip, (pdf(ip,iv), iv=1,nvar)
     ENDDO

     ip = nbins
     WRITE(21,1030) 1, 1, ip, (pdf(ip,iv), iv=1,nvar), (pdf(nbins+2,iv),iv=1,nvar)

1020 FORMAT(I5,2(1X,I5),L_FORMAT_MAX(1X,E12.3E3))
1030 FORMAT(I5,2(1X,I5),L_FORMAT_MAX(1X,E12.3E3),L_FORMAT_MAX(1X,E12.3E3))

     CLOSE(21)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE PDF3D_N


