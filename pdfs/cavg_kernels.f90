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
SUBROUTINE CAVG1V2D(ilim, nx,ny,nz, j, umin_ext,umax_ext,u, a, nbins,pdf,avg,  wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, j, nbins
  TREAL umin_ext,umax_ext
  TREAL, INTENT(IN)    :: u(nx,ny,nz), a(nx,ny,nz)
  TREAL, INTENT(OUT)   :: avg(nbins+2) ! Space at the end for min/max values in the sample variable
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

  ustep = (umax-umin) /M_REAL(nbins)    ! Calculate step in histogram
  avg(nbins+1) = umin +C_05_R *ustep    ! Calculate coordinate of histogram
  avg(nbins+2) = umax -C_05_R *ustep
  IF ( ustep .EQ. C_0_R ) ustep = C_1_R ! Just 1 point, prevent division by zero and force all in first bin

  ! -------------------------------------------------------------------
  ! Calculate Average and Histogram
  ! -------------------------------------------------------------------
  DO k = 1,nz
    DO i = 1,nx
      up = INT((u(i,j,k)-umin)/ustep) + 1
      IF ( ilim .EQ. 0 ) THEN
        IF ( up .LE. nbins .AND. up .GE. 1 ) THEN
          pdf(up) = pdf(up) + C_1_R
          avg(up) = avg(up) + a(i,j,k)
        ENDIF
      ELSE ! put last point in the last bin
        up = MIN(up,nbins)
        pdf(up) = pdf(up) + C_1_R
        avg(up) = avg(up) + a(i,j,k)
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
SUBROUTINE CAVG1V2D1G(ilim, nx,ny,nz, j, igate,gate, umin_ext,umax_ext,u, a, nbins,pdf,avg, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, j, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx,ny,nz), a(nx,ny,nz)
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

  ustep = (umax-umin) /M_REAL(nbins)    ! Calculate Step in Histogram
  avg(nbins+1) = umin +C_05_R *ustep    ! Calculate coordinate of histogram
  avg(nbins+2) = umax -C_05_R *ustep
  IF ( ustep .EQ. C_0_R ) ustep = C_1_R ! Just 1 point, prevent division by zero and force all in first bin

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  DO k = 1,nz
    DO i = 1,nx
      IF ( gate(i,j,k) .EQ. igate ) THEN
        up = INT((u(i,j,k)-umin)/ustep) + 1
        IF ( ilim .EQ. 0 ) THEN
          IF ( up .LE. nbins .AND. up .GE. 1 ) THEN
            pdf(up) = pdf(up) + C_1_R
            avg(up) = avg(up) + a(i,j,k)
          ENDIF
        ELSE ! put last point in the last bin
          up = MIN(up,nbins)
          pdf(up) = pdf(up) + C_1_R
          avg(up) = avg(up) + a(i,j,k)
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

  RETURN
END SUBROUTINE CAVG1V2D1G

!########################################################################
! Conditioned on 2 variables
!########################################################################
SUBROUTINE CAVG2V2D(nx,ny,nz, j, u,v, a, nbins,pdf,avg, wrk2d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, j, nbins(2)
  TREAL, INTENT(IN)    :: u(nx,ny,nz),v(nx,ny,nz), a(nx,ny,nz)
  TREAL, INTENT(OUT)   :: avg(nbins(1)*nbins(2) +2 +2*nbins(1)) ! Space at the end for min/max values of sample variable
  TREAL, INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))              ! nbins(2) should be greater than 2 for enough memory space
  TREAL, INTENT(INOUT) :: pdf(nbins(1)*nbins(2))

  ! -------------------------------------------------------------------
  TINTEGER i,k, up,vp, ip, offset
  TREAL umin,umax,ustep

#ifdef USE_MPI
  INTEGER ims_err, impi
  TREAL umin_p, umax_p
#endif

  ! ###################################################################
  pdf = C_0_R
  avg = C_0_R

  offset = nbins(1)*nbins(2) +2

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  ! First variable
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

  ustep = (umax-umin) /M_REAL(nbins(1))         ! Calculate step in histogram
  avg(nbins(1)*nbins(2)+1) = umin +C_05_R*ustep ! Calculate coordinate of histogram
  avg(nbins(1)*nbins(2)+2) = umax -C_05_R*ustep
  IF ( ustep .EQ. C_0_R) ustep = C_1_R          ! Just 1 point, prevent division by zero and force all in first bin

  ! Second variable
  ip = offset +1;    avg(ip) = v(1,j,1)
  ip = ip +nbins(1); avg(ip) = v(1,j,1)
  DO k = 1,nz
    DO i = 1,nx
      up = INT((u(i,j,k)-umin)/ustep) + 1
      up = MAX(1,MIN(up,nbins(1)))
      ip = offset +up;   avg(ip) = MIN(avg(ip),v(i,j,k))
      ip = ip +nbins(1); avg(ip) = MAX(avg(ip),v(i,j,k))
    ENDDO
  ENDDO
#ifdef USE_MPI
  impi = nbins(1)
  ip = offset +1;    CALL MPI_ALLREDUCE(avg(ip), wrk2d, impi, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
  avg(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
  ip = ip +nbins(1); CALL MPI_ALLREDUCE(avg(ip), wrk2d, impi, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  avg(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
#endif

  DO up = 1,nbins(1)                                        ! Calculate Step in Histogram
    wrk2d(up) = ( avg(offset+up+nbins(1)) -avg(offset+up) ) /M_REAL(nbins(2))
    IF ( wrk2d(up) .EQ. C_0_R ) wrk2d(up) = C_1_R           ! Just 1 point, prevent division by zero and force all in first bin
  ENDDO

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  DO k = 1,nz
    DO i = 1,nx
      up = INT((u(i,j,k)-umin)          /ustep    ) + 1
      up = MAX(1,MIN(up,nbins(1)))
      vp = INT((v(i,j,k)-avg(offset+up))/wrk2d(up)) + 1
      vp = MAX(1,MIN(vp,nbins(2)))
      ip = (vp-1)*nbins(1) +up
      pdf(ip) = pdf(ip) + C_1_R
      avg(ip) = avg(ip) + a(i,j,k)
    ENDDO
  ENDDO

  DO up = 1,nbins(1)                                        ! Calculate coordinate of histogram; I needed the minimum before
    ip = offset +up;   avg(ip) = avg(ip) +C_05_R*wrk2d(up)
    ip = ip +nbins(1); avg(ip) = avg(ip) -C_05_R*wrk2d(up)
  ENDDO

#ifdef USE_MPI
  impi = nbins(1)*nbins(2)
  CALL MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
  CALL MPI_ALLREDUCE(avg, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  avg(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
#endif

  DO ip = 1,nbins(1)*nbins(2)
    IF ( pdf(ip) .GT. C_0_R ) THEN ! Avg remains zero if there is no point in this interval
      avg(ip) = avg(ip) /pdf(ip)
    ENDIF
  ENDDO

  RETURN
END SUBROUTINE CAVG2V2D

!########################################################################
! Now the same, but using calculating the average over the whole array u
!########################################################################
SUBROUTINE CAVG2V3D(nx,ny,nz, u,v, a, nbins,pdf,avg, wrk2d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, nbins(2)
  TREAL, INTENT(IN)    :: u(nx*ny*nz),v(nx*ny*nz), a(nx*ny*nz)
  TREAL, INTENT(OUT)   :: avg(nbins(1)*nbins(2) +2 +2*nbins(1)) ! Space at the end for min/max values of sample variable
  TREAL, INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))              ! nbins(2) should be greater than 2 for enough memory space
  TREAL, INTENT(INOUT) :: pdf(nbins(1)*nbins(2))

  ! -------------------------------------------------------------------
  TINTEGER i, up,vp, ip, offset
  TREAL umin,umax,ustep
#ifdef USE_MPI
  INTEGER ims_err, impi
#endif

  ! ###################################################################
  pdf = C_0_R
  avg = C_0_R

  offset = nbins(1)*nbins(2) +2

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  ! First variable
  CALL MINMAX(nx,ny,nz, u, umin,umax)

  ustep = (umax-umin) /M_REAL(nbins(1))         ! Calculate step in histogram
  avg(nbins(1)*nbins(2)+1) = umin +C_05_R*ustep ! Calculate coordinate of histogram
  avg(nbins(1)*nbins(2)+2) = umax -C_05_R*ustep
  IF ( ustep .EQ. C_0_R) ustep = C_1_R          ! Just 1 point, prevent division by zero and force all in first bin

  ! Second variable
  ip = offset +1;    avg(ip) = v(1)
  ip = ip +nbins(1); avg(ip) = v(1)
  DO i = 1,nx*ny*nz
    up = INT((u(i)-umin)/ustep) + 1
    up = MAX(1,MIN(up,nbins(1)))
    ip = offset +up;   avg(ip) = MIN(avg(ip),v(i))
    ip = ip +nbins(1); avg(ip) = MAX(avg(ip),v(i))
  ENDDO
#ifdef USE_MPI
  impi = nbins(1)
  ip = offset +1;    CALL MPI_ALLREDUCE(avg(ip), wrk2d, impi, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
  avg(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
  ip = ip +nbins(1); CALL MPI_ALLREDUCE(avg(ip), wrk2d, impi, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  avg(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
#endif

  DO up = 1,nbins(1)                                        ! Calculate Step in Histogram
    wrk2d(up) = ( avg(offset+up+nbins(1)) -avg(offset+up) ) /M_REAL(nbins(2))
    IF ( wrk2d(up) .EQ. C_0_R ) wrk2d(up) = C_1_R           ! Just 1 point, prevent division by zero and force all in first bin
  ENDDO

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  DO i = 1,nx*ny*nz
    up = INT((u(i)-umin)          /ustep    ) + 1
    up = MAX(1,MIN(up,nbins(1)))
    vp = INT((v(i)-avg(offset+up))/wrk2d(up)) + 1
    vp = MAX(1,MIN(vp,nbins(2)))
    ip = (vp-1)*nbins(1) +up
    pdf(ip) = pdf(ip) + C_1_R
    avg(ip) = avg(ip) + a(i)
  ENDDO

  DO up = 1,nbins(1)                                        ! Calculate coordinate of histogram; I needed the minimum before
    ip = offset +up;   avg(ip) = avg(ip) +C_05_R*wrk2d(up)
    ip = ip +nbins(1); avg(ip) = avg(ip) -C_05_R*wrk2d(up)
  ENDDO

#ifdef USE_MPI
  impi = nbins(1)*nbins(2)
  CALL MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
  CALL MPI_ALLREDUCE(avg, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  avg(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
#endif

  DO ip = 1,nbins(1)*nbins(2)
    IF ( pdf(ip) .GT. C_0_R ) THEN ! Avg remains zero if there is no point in this interval
      avg(ip) = avg(ip) /pdf(ip)
    ENDIF
  ENDDO

  RETURN
END SUBROUTINE CAVG2V3D
