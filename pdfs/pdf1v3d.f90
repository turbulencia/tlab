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
SUBROUTINE PDF1V3D(ilim, imax,jmax,kmax, umin_ext, umax_ext, u, nbins, pdf, wrk1d)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, imax,jmax,kmax, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(imax*jmax*kmax)
  TREAL, INTENT(OUT)   :: pdf(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins)

! -------------------------------------------------------------------
  TINTEGER i, ip
  TREAL umin, umax, pdfstep
#ifdef USE_MPI
  INTEGER impi, ims_err
#endif

! ###################################################################
  pdf = C_0_R
  
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
  pdf(nbins+1) = umin + pdfstep/C_2_R
  pdf(nbins+2) = umax - pdfstep/C_2_R

! -------------------------------------------------------------------
! Calculate Histogram
! -------------------------------------------------------------------
  IF ( ABS(pdfstep) .GT. C_0_R ) THEN
     DO i = 1,imax*jmax*kmax
        ip = INT((u(i)-umin)/pdfstep) + 1
        IF ( ilim .EQ. 0 ) THEN
           IF ( ip .LE. nbins .AND. ip .GE. 1 ) THEN
              pdf(ip) = pdf(ip) + C_1_R
           ENDIF
        ELSE ! put last point in the last bin
           ip = MIN(ip,nbins)
           pdf(ip) = pdf(ip) + C_1_R
        ENDIF
     ENDDO

#ifdef USE_MPI
     impi = nbins
     CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
     pdf(1:nbins) = wrk1d(1:nbins)
#endif

  ENDIF

  RETURN
END SUBROUTINE PDF1V3D

