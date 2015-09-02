!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/07/11 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Computes the conditional PDF P(v|u) using the whole 3D array. It simply
!# calculates the PDF of v for points around ucenter.
!#
!########################################################################
!# ARGUMENTS 
!#
!# ilim     In    If set =0, externally forced through vmin_ext/vmax_ext
!#                If not, calculate locally the min/max
!#
!########################################################################
SUBROUTINE CPDF2V3D(inorm, ilim, imax, jmax, kmax, vmin_ext, vmax_ext,&
     u, v, ucenter, udelta, nbins, y, pdf, wrk1d, nsample)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER inorm, ilim
  TINTEGER imax, jmax, kmax
  TREAL vmin_ext, vmax_ext, ucenter, udelta
  TINTEGER nbins, nsample
  TREAL u(*)
  TREAL v(*)
  TREAL y(nbins)
  TREAL pdf(nbins)
  TREAL wrk1d(nbins)

! -------------------------------------------------------------------
  TINTEGER jp, i
  TREAL umin, umax, vmin, vmax
  TREAL pdfstep, pnorm
#ifdef USE_MPI
  TINTEGER nsample_mpi
  INTEGER impi, ims_err
  TREAL vmin_p, vmax_p
#endif

! ###################################################################
  nsample = 0
  DO jp = 1, nbins
     pdf(jp) = C_0_R
  ENDDO

! ###################################################################
  IF ( ABS(udelta) .GT. C_0_R ) THEN
     umin = ucenter - udelta*C_05_R
     umax = ucenter + udelta*C_05_R

! -------------------------------------------------------------------
! Calculate the conditional min/max values
! -------------------------------------------------------------------
     IF ( ilim .EQ. 0 ) THEN
        vmin = vmin_ext
        vmax = vmax_ext

     ELSE
        vmin = C_BIG_R
        vmax =-C_BIG_R
        DO i = 1, imax*jmax*kmax
           IF ( u(i) .GE. umin .AND. u(i) .LE. umax ) THEN
              vmin = MIN(vmin, v(i))
              vmax = MAX(vmax, v(i))
           ENDIF
        ENDDO

#ifdef USE_MPI
        CALL MPI_ALLREDUCE(vmin, vmin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
        vmin = vmin_p
        CALL MPI_ALLREDUCE(vmax, vmax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        vmax = vmax_p
#endif

     ENDIF

! Calculate Step in Histogram
     pdfstep = (vmax-vmin)/M_REAL(nbins)

! Calculate x coordinate of histogram
     y(1) = vmin + pdfstep*C_05_R
     DO jp=2, nbins
        y(jp) = y(jp-1) + pdfstep
     ENDDO

! -------------------------------------------------------------------
! Calculate Histogram
! -------------------------------------------------------------------
     IF ( ABS(pdfstep) .GT. C_0_R ) THEN
        DO i = 1, imax*jmax*kmax
           IF ( u(i) .GE. umin .AND. u(i) .LE. umax ) THEN
              jp = INT((v(i)-vmin)/pdfstep) + 1
              IF ( ilim .EQ. 0 ) THEN
                 IF ( jp .LE. nbins .AND. jp .GE. 1 ) THEN
                    pdf(jp) = pdf(jp) + C_1_R
                    nsample = nsample + 1
                 ENDIF
              ELSE ! put last point in the last bin
                 jp = MIN(jp,nbins)
                 pdf(jp) = pdf(jp) + C_1_R
                 nsample = nsample + 1
              ENDIF
           ENDIF
        ENDDO

#ifdef USE_MPI
! Sum all number of points
        nsample_mpi = nsample
        CALL MPI_ALLREDUCE(nsample_mpi,nsample, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)

        impi = nbins
        CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        DO jp=1, nbins
           pdf(jp) = wrk1d(jp)
        ENDDO
#endif
! normalize
        IF ( inorm .EQ. 1 ) THEN
           IF ( nsample .GT. 0 ) THEN
              pnorm = C_1_R/(M_REAL(nsample)*pdfstep)
              DO jp=1, nbins
                 pdf(jp) = pdf(jp)*pnorm
              ENDDO
           ELSE
              DO jp=1, nbins
                 pdf(jp) = C_0_R
                 y(jp)   = C_0_R
              ENDDO
           ENDIF
        ENDIF

     ENDIF
  ENDIF

  RETURN
END SUBROUTINE CPDF2V3D

