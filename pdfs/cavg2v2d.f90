!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/07/11 - J.P. Mellado
!#              Cleaned. 
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate conditional average <V|U> , 2 variables and along planes.
!# Array nsample2 is used in PARALLEL mode.
!# Same as CAVG3V2D but w/o gate field.
!#
!########################################################################
!# ARGUMENTS 
!#
!# imom     In    Moment order
!# ilim     In    If set =0, externally forced through umin_ext/umax_ext
!#                If not, calculate locally the min/max of conditioned u
!#
!########################################################################
SUBROUTINE CAVG2V2D&
     (ilim, imax, jmax, kmax, j, umin_ext, umax_ext, u, v, nbins, imom, x, cavg, nsample, wrk1d)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim
  TINTEGER imax, jmax, kmax, j, nbins, imom
  TREAL u(imax,jmax,kmax), umin_ext, umax_ext
  TREAL v(imax,jmax,kmax)
  TREAL cavg(nbins), wrk1d(nbins)
  TREAL x(nbins)
  TINTEGER nsample(nbins,2)

! -------------------------------------------------------------------
  TINTEGER ip, i, k
  TREAL umin, umax
  TREAL pdfstep

#ifdef USE_MPI
  INTEGER ims_err, impi
  TREAL umin_p, umax_p
#endif

! ###################################################################
  DO ip = 1, nbins
     nsample(ip,1) = 0
     cavg(ip)      = C_0_R
     x(ip)         = C_0_R
  ENDDO

! -------------------------------------------------------------------
! Calculate Minimum and Maximum
! -------------------------------------------------------------------
  IF ( ilim .EQ. 0 ) THEN
     umin = umin_ext
     umax = umax_ext

  ELSE
     umin = C_BIG_R
     umax =-C_BIG_R
     DO k = 1,kmax
        DO i = 1,imax
           umin = MIN(umin, u(i,j,k))
           umax = MAX(umax, u(i,j,k))
        ENDDO
     ENDDO

#ifdef USE_MPI
     CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
     umin = umin_p
     CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
     umax = umax_p
#endif
     
  ENDIF

! Calculate Step in Histogram
  pdfstep = (umax-umin)/M_REAL(nbins)

! -------------------------------------------------------------------
! Calculate Conditional Average
! -------------------------------------------------------------------
  IF ( pdfstep .GT. C_0_R ) THEN
     DO k = 1,kmax
        DO i = 1,imax
           ip = INT((u(i,j,k)-umin)/pdfstep) + 1
           IF ( ilim .EQ. 0 ) THEN
              IF ( ip .LE. nbins .AND. ip .GE. 1 ) THEN
                 cavg(ip)      = cavg(ip)      + v(i,j,k)**imom
                 nsample(ip,1) = nsample(ip,1) + 1
              ENDIF
           ELSE ! put last point in the last bin
              ip = MIN(ip,nbins)
              cavg(ip)      = cavg(ip)      + v(i,j,k)**imom
              nsample(ip,1) = nsample(ip,1) + 1
           ENDIF

        ENDDO
     ENDDO

#ifdef USE_MPI
! Sum all number of points
     impi = nbins
     CALL MPI_ALLREDUCE(cavg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
     CALL MPI_ALLREDUCE(nsample(1,1), nsample(1,2), impi, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)

     DO ip = 1, nbins
        cavg(ip)      = wrk1d(ip)
        nsample(ip,1) = nsample(ip,2)
     ENDDO
#endif

     DO ip = 1, nbins
        x(ip) = umin + pdfstep*(2*ip-1)/C_2_R
        IF ( nsample(ip,1) .EQ. 0 ) THEN
           cavg(ip) = C_0_R
        ELSE
           cavg(ip) = cavg(ip)/M_REAL(nsample(ip,1))
        ENDIF
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE CAVG2V2D

