#include "types.h"

MODULE PDFS
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: PDF1V2D, PDF1V2D1G, PDF2V2D, PDF_ANALIZE

  TINTEGER i,k, up,vp, ip, offset
  TREAL umin,umax,ustep

#ifdef USE_MPI
#include "mpif.h"
  INTEGER ims_err, impi
  TREAL umin_p, umax_p
#endif

CONTAINS
  !########################################################################
  !#
  !# Calculate the PDF over plane of an array u using nbins bins.
  !#
  !# ilim     In    0, externally forced through umin_ext/umax_ext
  !#                otherwise, calculate locally the min/max
  !#
  !########################################################################
  SUBROUTINE PDF1V2D(ilim, nx,ny,nz, j, umin_ext,umax_ext,u, nbins,pdf, wrk1d, a,avg)
    IMPLICIT NONE

    TINTEGER ilim, nx,ny,nz, j, nbins
    TREAL umin_ext,umax_ext
    TREAL, INTENT(IN   ) :: u(nx,ny,nz)
    TREAL, INTENT(  OUT) :: pdf(nbins+2)            ! Space at the end for min/max bins of u
    TREAL, INTENT(INOUT) :: wrk1d(nbins)
    TREAL, OPTIONAL      :: a(nx,ny,nz), avg(nbins) ! For conditional average, if needed

    ! ###################################################################
    pdf = C_0_R
    IF ( PRESENT(avg) ) avg = C_0_R

    ! -------------------------------------------------------------------
    ! Calculate Minimum and Maximum
    ! -------------------------------------------------------------------
    IF ( ilim == 0 ) THEN
      umin = umin_ext
      umax = umax_ext

    ELSE
      umin = u(1,j,1)
      umax = u(1,j,1)
      DO k = 1,nz
        DO i = 1,nx
          umin = MIN(umin, u(i,j,k))
          umax = MAX(umax, u(i,j,k))
        END DO
      END DO

#ifdef USE_MPI
      CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
      umin = umin_p
      CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
      umax = umax_p
#endif

    END IF

    ustep = (umax-umin) /M_REAL(nbins)    ! Calculate step in histogram
    pdf(nbins+1) = umin +C_05_R *ustep    ! Calculate coordinate of histogram
    pdf(nbins+2) = umax -C_05_R *ustep
    IF ( ustep == C_0_R ) ustep = C_1_R   ! Just 1 point, prevent division by zero and force all in first bin

    ! -------------------------------------------------------------------
    ! Calculate Histogram
    ! -------------------------------------------------------------------
    DO k = 1,nz
      DO i = 1,nx
        up = INT((u(i,j,k)-umin)/ustep) + 1
        IF ( ilim == 0 ) THEN
          IF ( up <= nbins .AND. up >= 1 ) THEN
            pdf(up) = pdf(up) + C_1_R
            IF ( PRESENT(a) .AND. PRESENT(avg) ) avg(up) = avg(up) + a(i,j,k)
          END IF
        ELSE ! put last point in the last bin
          up = MIN(up,nbins)
          pdf(up) = pdf(up) + C_1_R
          IF ( PRESENT(a) .AND. PRESENT(avg) ) avg(up) = avg(up) + a(i,j,k)
        END IF
      END DO
    END DO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
    IF ( PRESENT(avg) ) THEN
      CALL MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
      avg(1:nbins) = wrk1d(1:nbins)
    END IF
#endif

    IF ( PRESENT(avg) ) THEN              ! Save avg data in pdf array
      DO up = 1,nbins
        IF ( pdf(up) > C_0_R ) THEN       ! Avg remains zero if there is no point in this interval
          avg(up) = avg(up) /pdf(up)
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE PDF1V2D

  !########################################################################
  !#
  !# Conditioned on the intermittency field gate. Same as before, but with a
  !# conditional inside the loop.
  !#
  !# igate  In   Level of the gate signal to use as intermittency function
  !#
  !########################################################################
  SUBROUTINE PDF1V2D1G(ilim, nx,ny,nz, j, igate,gate, umin_ext,umax_ext,u, nbins,pdf, wrk1d, a,avg)
    IMPLICIT NONE

    TINTEGER ilim, nx,ny,nz, j, nbins
    TREAL umin_ext, umax_ext
    TREAL, INTENT(IN   ) :: u(nx,ny,nz)
    TREAL, INTENT(  OUT) :: pdf(nbins+2)            ! Space at the end for min/max bins of u
    TREAL, INTENT(INOUT) :: wrk1d(nbins)
    INTEGER(1),INTENT(IN):: gate(nx,ny,nz), igate
    TREAL, OPTIONAL      :: a(nx,ny,nz), avg(nbins) ! For conditional average, if needed

    ! ###################################################################
    pdf = C_0_R
    IF ( PRESENT(avg) ) avg = C_0_R

    ! -------------------------------------------------------------------
    ! Calculate Minimum and Maximum
    ! -------------------------------------------------------------------
    IF ( ilim == 0 ) THEN
      umin = umin_ext
      umax = umax_ext

    ELSE
      umin = u(1,j,1)
      umax = u(1,j,1)
      DO k = 1,nz
        DO i = 1,nx
          IF ( gate(i,j,k) == igate ) THEN
            umin = MIN(umin, u(i,j,k))
            umax = MAX(umax, u(i,j,k))
          END IF
        END DO
      END DO

#ifdef USE_MPI
      CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
      umin = umin_p
      CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
      umax = umax_p
#endif

    END IF

    ustep = (umax-umin) /M_REAL(nbins)    ! Calculate step in histogram
    pdf(nbins+1) = umin +C_05_R *ustep    ! Calculate coordinate of histogram
    pdf(nbins+2) = umax -C_05_R *ustep
    IF ( ustep == C_0_R ) ustep = C_1_R   ! Just 1 point, prevent division by zero and force all in first bin

    ! -------------------------------------------------------------------
    ! Calculate Histogram
    ! -------------------------------------------------------------------
    DO k = 1,nz
      DO i = 1,nx
        IF ( gate(i,j,k) == igate ) THEN
          up = INT((u(i,j,k)-umin)/ustep) + 1
          IF ( ilim == 0 ) THEN
            IF ( up <= nbins .AND. up >= 1 ) THEN
              pdf(up) = pdf(up) + C_1_R
              IF ( PRESENT(a) .AND. PRESENT(avg) ) avg(up) = avg(up) + a(i,j,k)
            END IF
          ELSE ! put last point in the last bin
            up = MIN(up,nbins)
            pdf(up) = pdf(up) + C_1_R
            IF ( PRESENT(a) .AND. PRESENT(avg) ) avg(up) = avg(up) + a(i,j,k)
          END IF
        END IF

      END DO
    END DO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
    IF ( PRESENT(avg) ) THEN
      CALL MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
      avg(1:nbins) = wrk1d(1:nbins)
    END IF
#endif

    IF ( PRESENT(avg) ) THEN              ! Save avg data in pdf array
      DO up = 1,nbins
        IF ( pdf(up) > C_0_R ) THEN       ! Avg remains zero if there is no point in this interval
          avg(up) = avg(up) /pdf(up)
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE PDF1V2D1G

  !########################################################################
  ! Joint PDFs
  !########################################################################
  SUBROUTINE PDF2V2D(nx,ny,nz, j, u,v, nbins,pdf, wrk2d, a,avg)
    IMPLICIT NONE

    TINTEGER nx,ny,nz, j, nbins(2)
    TREAL, INTENT(IN   ) :: u(nx,ny,nz), v(nx,ny,nz)
    TREAL, INTENT(  OUT) :: pdf(nbins(1)*nbins(2) +2 +2*nbins(1)) ! Space at the end for min/max bins of u,v
    TREAL, INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))              ! nbins(2) should be greater than 2 for enough memory space
    TREAL, OPTIONAL      :: a(nx,ny,nz), avg(nbins(1)*nbins(2))   ! For conditional average, if needed

    ! ###################################################################
    pdf = C_0_R
    IF ( PRESENT(avg) ) avg = C_0_R

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
      END DO
    END DO

#ifdef USE_MPI
    CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
    umin = umin_p
    CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
    umax = umax_p
#endif

    ustep = (umax-umin) /M_REAL(nbins(1))         ! Calculate step in histogram
    pdf(nbins(1)*nbins(2)+1) = umin +C_05_R*ustep ! Calculate coordinate of histogram
    pdf(nbins(1)*nbins(2)+2) = umax -C_05_R*ustep
    IF ( ustep == C_0_R) ustep = C_1_R            ! Just 1 point, prevent division by zero and force all in first bin

    ! Second variable
    ip = offset +1;    pdf(ip) = v(1,j,1)
    ip = ip +nbins(1); pdf(ip) = v(1,j,1)
    DO k = 1,nz
      DO i = 1,nx
        up = INT((u(i,j,k)-umin)/ustep) + 1
        up = MAX(1,MIN(up,nbins(1)))
        ip = offset +up;   pdf(ip) = MIN(pdf(ip),v(i,j,k))
        ip = ip +nbins(1); pdf(ip) = MAX(pdf(ip),v(i,j,k))
      END DO
    END DO
#ifdef USE_MPI
    impi = nbins(1)
    ip = offset +1;    CALL MPI_ALLREDUCE(pdf(ip), wrk2d, impi, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
    pdf(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
    ip = ip +nbins(1); CALL MPI_ALLREDUCE(pdf(ip), wrk2d, impi, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
    pdf(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
#endif

    DO up = 1,nbins(1)                                          ! Calculate step in histogram
      wrk2d(up) = ( pdf(offset+up+nbins(1)) -pdf(offset+up) ) /M_REAL(nbins(2))
      IF ( wrk2d(up) == C_0_R ) wrk2d(up) = C_1_R               ! Just 1 point, prevent division by zero and force all in first bin
    END DO

    ! -------------------------------------------------------------------
    ! Calculate Histogram
    ! -------------------------------------------------------------------
    DO k = 1,nz
      DO i = 1,nx
        up = INT((u(i,j,k)-umin)          /ustep    ) + 1
        up = MAX(1,MIN(up,nbins(1)))
        vp = INT((v(i,j,k)-pdf(offset+up))/wrk2d(up)) + 1
        vp = MAX(1,MIN(vp,nbins(2)))
        ip = (vp-1)*nbins(1) +up
        pdf(ip) = pdf(ip) + C_1_R
        IF ( PRESENT(a) .AND. PRESENT(avg) ) avg(ip) = avg(ip) + a(i,j,k)
      END DO
    END DO

    DO up = 1,nbins(1)                                        ! Calculate coordinate of histogram; I needed the minimum before
      ip = offset +up;   pdf(ip) = pdf(ip) +C_05_R*wrk2d(up)
      ip = ip +nbins(1); pdf(ip) = pdf(ip) -C_05_R*wrk2d(up)
    END DO

#ifdef USE_MPI
    impi = nbins(1)*nbins(2)
    CALL MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
    IF ( PRESENT(avg) ) THEN
      CALL MPI_ALLREDUCE(avg, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
      avg(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
    END IF
#endif

    IF ( PRESENT(avg) ) THEN              ! Save avg data in pdf array
      DO ip = 1,nbins(1)*nbins(2)
        IF ( pdf(ip) > C_0_R ) THEN       ! Avg remains zero if there is no point in this interval
          avg(ip) = avg(ip) /pdf(ip)
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE PDF2V2D

  !########################################################################
  !# Recalculating max/min for a given relative threshold plim.
  !# Adding nplim, the count of points above threshold.
  !# BCs flag ibc to drop extreme points.
  !########################################################################
  SUBROUTINE PDF_ANALIZE(ibc, nbins,pdf, umin_ext,umax_ext, plim,nplim)
    IMPLICIT NONE

    TINTEGER, INTENT(IN   ) :: ibc,nbins
    TREAL,    INTENT(IN   ) :: pdf(nbins+2), plim
    TREAL,    INTENT(INOUT) :: umin_ext,umax_ext
    TINTEGER, INTENT(  OUT) :: nplim

    ! -------------------------------------------------------------------
    TINTEGER upmin,upmax
    TREAL pdf_threshold

    ! ###################################################################
    ustep = ( pdf(nbins+2) -pdf(nbins+1) ) /M_REAL(nbins-1)

    upmin = 1                                       ! eliminate the outer bins according to BCs
    upmax = nbins
    IF ( ibc == 1 .OR. ibc == 3 ) THEN
      upmin = upmin + 1
    END IF
    IF ( ibc == 2 .OR. ibc == 3 ) THEN
      upmax = upmax - 1
    END IF
    umin_ext = pdf(nbins+1) -C_05_R*ustep +ustep*M_REAL(upmin-1)
    umax_ext = pdf(nbins+1) -C_05_R*ustep +ustep*M_REAL(upmax)

    pdf_threshold = plim *MAXVAL(pdf(upmin:upmax))  ! get absolute threshold

    nplim = 0                                       ! count the number of points above threshold
    DO up = upmin,upmax
      IF ( pdf(up) > pdf_threshold ) THEN
        nplim = nplim + 1
      END IF
    END DO

    DO up = upmin,upmax                             ! elmininate smallest u-values if their probability is below threshold
      IF ( pdf(up) > pdf_threshold ) THEN
        umin_ext = pdf(nbins+1) -C_05_R*ustep +ustep*M_REAL(up-1)
        EXIT
      END IF
    END DO

    DO up = upmax,upmin,-1                          ! elmininate largest u-values if their probability is below threshold
      IF ( pdf(up) > pdf_threshold ) THEN
        umax_ext = pdf(nbins+1) -C_05_R*ustep +ustep*M_REAL(up)
        EXIT
      END IF
    END DO

    RETURN
  END SUBROUTINE PDF_ANALIZE

END MODULE PDFS
