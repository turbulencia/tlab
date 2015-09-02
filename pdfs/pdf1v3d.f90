!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/07/11 - J.P. Mellado
!#              Cleaned. Similar to PDF1V2D
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the PDF over an array u using nbins bins.
!#
!########################################################################
!# ARGUMENTS 
!#
!# ilim     In    If set =0, externally forced through umin_ext/umax_ext
!#                If not, calculate locally the min/max
!#
!########################################################################
SUBROUTINE PDF1V3D(inorm, ilim, imax, jmax, kmax, &
     umin_ext, umax_ext, u, nbins, xpdfmin, xpdfmax, pdf, wrk1d)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER inorm, ilim
  TINTEGER imax, jmax, kmax
  TREAL umin_ext, umax_ext
  TINTEGER nbins
  TREAL u(*)
  TREAL xpdfmin, xpdfmax
  TREAL pdf(nbins)
  TREAL wrk1d(nbins)

! -------------------------------------------------------------------
  TINTEGER i
  TREAL umin, umax
  TREAL pdfstep, pnorm
  TINTEGER ip, nsample
#ifdef USE_MPI
  TINTEGER nsample_mpi
  INTEGER impi, ims_err
  TREAL umin_p, umax_p
#endif

! ###################################################################
  nsample = 0
  DO ip = 1,nbins
     pdf(ip) = C_0_R
  ENDDO

! -------------------------------------------------------------------
! Calculate Minimum and Maximum
! -------------------------------------------------------------------
  IF ( ilim .EQ. 0 ) THEN
     umin = umin_ext
     umax = umax_ext

  ELSE
     CALL MINMAX(imax,jmax,kmax, u, umin,umax)

  ENDIF

! Calculate Step in Histogram
  pdfstep = (umax-umin)/M_REAL(nbins)

! Calculate x coordinate of histogram
  xpdfmin = umin + pdfstep/C_2_R
  xpdfmax = umin + pdfstep*(2*nbins-1)/C_2_R

! -------------------------------------------------------------------
! Calculate Histogram
! -------------------------------------------------------------------
  IF ( ABS(pdfstep) .GT. C_0_R ) THEN
     DO i = 1,imax*jmax*kmax
        ip = INT((u(i)-umin)/pdfstep) + 1
        IF ( ilim .EQ. 0 ) THEN
           IF ( ip .LE. nbins .AND. ip .GE. 1 ) THEN
              pdf(ip) = pdf(ip) + C_1_R
              nsample = nsample + 1
           ENDIF
        ELSE ! put last point in the last bin
           ip = MIN(ip,nbins)
           pdf(ip) = pdf(ip) + C_1_R
           nsample = nsample + 1
        ENDIF
     ENDDO

#ifdef USE_MPI
! Sum all number of points
     nsample_mpi = nsample
     CALL MPI_ALLREDUCE(nsample_mpi, nsample, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)

     impi = nbins
     CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
     DO ip = 1,nbins
        pdf(ip) = wrk1d(ip)
     ENDDO
#endif
! normalize
     IF ( inorm .EQ. 1 ) THEN
        IF ( nsample .GT. 0 ) THEN
           pnorm = C_1_R/(M_REAL(nsample)*pdfstep)
           DO ip = 1,nbins
              pdf(ip) = pdf(ip)*pnorm
           ENDDO
        ELSE
           DO ip = 1,nbins
              pdf(ip) = C_0_R
              xpdfmin = C_0_R
              xpdfmax = C_0_R
           ENDDO
        ENDIF
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE PDF1V3D

