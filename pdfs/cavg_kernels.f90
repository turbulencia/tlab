#include "types.h"

!########################################################################
!#
!# Calculate the average of v over plane conditioned on an array u using nbins bins.
!# This routine follows pdf1v2d, ...
!#
!# ilim     In    0, externally forced through umin_ext/umax_ext
!#                otherwise, calculate locally the min/max
!#
!########################################################################
SUBROUTINE CAVG1V2D(ilim, nx,ny,nz, j, umin_ext,umax_ext, u,v, nbins,pdf,avg,  wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, j, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx,ny,nz), v(nx,ny,nz)
  TREAL, INTENT(OUT)   :: avg(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins), pdf(nbins)

  ! -------------------------------------------------------------------
  TINTEGER i,k, up
  TREAL umin,umax,ustep

#ifdef USE_MPI
  INTEGER ims_err, impi
  TREAL umin_p, umax_p
#endif

  ! ###################################################################
  pdf = C_0_R
  avg = C_0_R

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  IF ( ilim .EQ. 0 ) THEN
    umin = umin_ext
    umax = umax_ext

  ELSE
    umin = u(1,j,1)
    umax = u(1,j,1)
    DO k = 1,nz
      DO i = 1,nx
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
  ustep = (umax-umin)/M_REAL(nbins)

  ! Calculate x coordinate of histogram
  avg(nbins+1) = umin + ustep/C_2_R
  avg(nbins+2) = umax - ustep/C_2_R

  ! -------------------------------------------------------------------
  ! Calculate Average and Histogram
  ! -------------------------------------------------------------------
  IF ( ABS(ustep) .GT. C_0_R ) THEN
    DO k = 1,nz
      DO i = 1,nx
        up = INT((u(i,j,k)-umin)/ustep) + 1
        IF ( ilim .EQ. 0 ) THEN
          IF ( up .LE. nbins .AND. up .GE. 1 ) THEN
            pdf(up) = pdf(up) + C_1_R
            avg(up) = avg(up) + v(i,j,k)
          ENDIF
        ELSE ! put last point in the last bin
          up = MIN(up,nbins)
          pdf(up) = pdf(up) + C_1_R
          avg(up) = avg(up) + v(i,j,k)
        ENDIF
      ENDDO
    ENDDO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
    CALL MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    avg(1:nbins) = wrk1d(1:nbins)
#endif

    DO up = 1,nbins
      IF ( pdf(up) .GT. C_0_R ) THEN ! Avg remains zero if there is no point in this interval
        avg(up) = avg(up) /pdf(up)
      ENDIF
    ENDDO

  ENDIF

  RETURN
END SUBROUTINE CAVG1V2D

!########################################################################
!#
!# Conditioned on the intermittency field gate. Same as before, but with a
!# conditional inside the loop.
!#
!# igate  In   Level of the gate signal to use as intermittency function
!#
!########################################################################
SUBROUTINE CAVG1V2D1G(ilim, nx,ny,nz, j, igate,gate, umin_ext,umax_ext,u, v, nbins,pdf,avg, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, j, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx,ny,nz), v(nx,ny,nz)
  TREAL, INTENT(OUT)   :: avg(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins), pdf(nbins)
  INTEGER(1),INTENT(IN):: gate(nx,ny,nz), igate

  ! -------------------------------------------------------------------
  TINTEGER i, k, up
  TREAL umin, umax, ustep

#ifdef USE_MPI
  INTEGER ims_err, impi
  TREAL umin_p, umax_p
#endif

  ! ###################################################################
  pdf = C_0_R
  avg = C_0_R

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  IF ( ilim .EQ. 0 ) THEN
    umin = umin_ext
    umax = umax_ext

  ELSE
    umin = u(1,j,1)
    umax = u(1,j,1)
    DO k = 1,nz
      DO i = 1,nx
        IF ( gate(i,j,k) .EQ. igate ) THEN
          umin = MIN(umin, u(i,j,k))
          umax = MAX(umax, u(i,j,k))
        ENDIF
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
  ustep = (umax-umin)/M_REAL(nbins)

  ! Calculate x coordinate of histogram
  pdf(nbins+1) = umin + ustep/C_2_R
  pdf(nbins+2) = umax - ustep/C_2_R

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  IF ( ABS(ustep) .GT. C_0_R ) THEN
    DO k = 1,nz
      DO i = 1,nx
        IF ( gate(i,j,k) .EQ. igate ) THEN
          up = INT((u(i,j,k)-umin)/ustep) + 1
          IF ( ilim .EQ. 0 ) THEN
            IF ( up .LE. nbins .AND. up .GE. 1 ) THEN
              pdf(up) = pdf(up) + C_1_R
              avg(up) = avg(up) + v(i,j,k)
            ENDIF
          ELSE ! put last point in the last bin
            up = MIN(up,nbins)
            pdf(up) = pdf(up) + C_1_R
            avg(up) = avg(up) + v(i,j,k)
          ENDIF
        ENDIF

      ENDDO
    ENDDO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
    CALL MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    avg(1:nbins) = wrk1d(1:nbins)
#endif

    DO up = 1,nbins
      IF ( pdf(up) .GT. C_0_R ) THEN ! Avg remains zero if there is no point in this interval
        avg(up) = avg(up) /pdf(up)
      ENDIF
    ENDDO

  ENDIF

  RETURN
END SUBROUTINE CAVG1V2D1G

!########################################################################
! Now the same, but using calculating the PDF over the whole array u
!########################################################################
SUBROUTINE CAVG1V3D(ilim, nx,ny,nz, umin_ext,umax_ext,u, v, nbins,pdf,avg, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx*ny*nz), v(nx*ny*nz)
  TREAL, INTENT(OUT)   :: avg(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins), pdf(nbins)

  ! -------------------------------------------------------------------
  TINTEGER i, up
  TREAL umin, umax, ustep
#ifdef USE_MPI
  INTEGER impi, ims_err
#endif

  ! ###################################################################
  pdf = C_0_R
  avg = C_0_R

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  IF ( ilim .EQ. 0 ) THEN
    umin = umin_ext
    umax = umax_ext

  ELSE
    CALL MINMAX(nx,ny,nz, u, umin,umax)

  ENDIF

  ! Calculate Step in Histogram
  ustep = (umax-umin)/M_REAL(nbins)

  ! Calculate x coordinate of histogram
  pdf(nbins+1) = umin + ustep/C_2_R
  pdf(nbins+2) = umax - ustep/C_2_R

  ! -------------------------------------------------------------------
! Calculate Histogram
  ! -------------------------------------------------------------------
  IF ( ABS(ustep) .GT. C_0_R ) THEN
    DO i = 1,nx*ny*nz
      up = INT((u(i)-umin)/ustep) + 1
      IF ( ilim .EQ. 0 ) THEN
        IF ( up .LE. nbins .AND. up .GE. 1 ) THEN
          pdf(up) = pdf(up) + C_1_R
          avg(up) = avg(up) + v(i)
        ENDIF
      ELSE ! put last point in the last bin
        up = MIN(up,nbins)
        pdf(up) = pdf(up) + C_1_R
        avg(up) = avg(up) + v(i)
      ENDIF
    ENDDO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
    CALL MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    avg(1:nbins) = wrk1d(1:nbins)
#endif

    DO up = 1,nbins
      IF ( pdf(up) .GT. C_0_R ) THEN ! Avg remains zero if there is no point in this interval
        avg(up) = avg(up) /pdf(up)
      ENDIF
    ENDDO

  ENDIF

  RETURN
END SUBROUTINE CAVG1V3D

!########################################################################
! And the conditioned one
!########################################################################
SUBROUTINE CAVG1V3D1G(ilim, nx,ny,nz, igate,gate, umin_ext,umax_ext,u, v, nbins,pdf,avg, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx*ny*nz), v(nx*ny*nz)
  TREAL, INTENT(OUT)   :: avg(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins), pdf(nbins)
  INTEGER(1),INTENT(IN):: gate(nx*ny*nz), igate

  ! -------------------------------------------------------------------
  TINTEGER i, up
  TREAL umin, umax, ustep

#ifdef USE_MPI
  INTEGER impi, ims_err
  TREAL umin_p, umax_p
#endif

  ! ###################################################################
  pdf = C_0_R
  avg = C_0_R

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  IF ( ilim .EQ. 0 ) THEN
    umin = umin_ext
    umax = umax_ext

  ELSE
    umin = u(1)
    umax = u(1)
    DO i = 1,nx*ny*nz
      IF ( gate(i) .EQ. igate ) THEN
        umin = MIN(umin, u(i))
        umax = MAX(umax, u(i))
      ENDIF
    ENDDO

#ifdef USE_MPI
    CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
    umin = umin_p
    CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
    umax = umax_p
#endif

  ENDIF

  ! Calculate Step in Histogram
  ustep = (umax-umin)/M_REAL(nbins)

  ! Calculate x coordinate of histogram
  pdf(nbins+1) = umin + ustep/C_2_R
  pdf(nbins+2) = umax - ustep/C_2_R

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  IF ( ABS(ustep) .GT. C_0_R ) THEN
    DO i = 1,nx*ny*nz
      IF ( gate(i) .EQ. igate ) THEN
        up = INT((u(i)-umin)/ustep) + 1
        IF ( ilim .EQ. 0 ) THEN
          IF ( up .LE. nbins .AND. up .GE. 1 ) THEN
            pdf(up) = pdf(up) + C_1_R
            avg(up) = avg(up) + v(i)
          ENDIF
        ELSE ! put last point in the last bin
          up = MIN(up,nbins)
          pdf(up) = pdf(up) + C_1_R
          avg(up) = avg(up) + v(i)
        ENDIF
      ENDIF

    ENDDO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
    CALL MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    avg(1:nbins) = wrk1d(1:nbins)
#endif

    DO up = 1,nbins
      IF ( pdf(up) .GT. C_0_R ) THEN ! Avg remains zero if there is no point in this interval
        avg(up) = avg(up) /pdf(up)
      ENDIF
    ENDDO

  ENDIF

  RETURN
END SUBROUTINE CAVG1V3D1G
