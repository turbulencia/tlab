#include "types.h"

!########################################################################
!#
!# Calculate the PDF over plane of an array u using nbins bins.
!#
!# ilim     In    0, externally forced through umin_ext/umax_ext
!#                otherwise, calculate locally the min/max
!#
!########################################################################
SUBROUTINE PDF1V2D(ilim, nx,ny,nz, j, umin_ext,umax_ext,u, nbins,pdf, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, j, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx,ny,nz)
  TREAL, INTENT(OUT)   :: pdf(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins)

  ! -------------------------------------------------------------------
  TINTEGER i,k, up
  TREAL umin,umax,ustep

#ifdef USE_MPI
  INTEGER ims_err, impi
  TREAL umin_p, umax_p
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
  pdf(nbins+1) = umin + ustep/C_2_R
  pdf(nbins+2) = umax - ustep/C_2_R

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  IF ( ABS(ustep) .GT. C_0_R ) THEN
    DO k = 1,nz
      DO i = 1,nx
        up = INT((u(i,j,k)-umin)/ustep) + 1
        IF ( ilim .EQ. 0 ) THEN
          IF ( up .LE. nbins .AND. up .GE. 1 ) THEN
            pdf(up) = pdf(up) + C_1_R
          ENDIF
        ELSE ! put last point in the last bin
          up = MIN(up,nbins)
          pdf(up) = pdf(up) + C_1_R
        ENDIF
      ENDDO
    ENDDO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
#endif

  ENDIF

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
SUBROUTINE PDF1V2D1G(ilim, nx,ny,nz, j, igate,gate, umin_ext,umax_ext,u, nbins,pdf, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, j, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx,ny,nz)
  TREAL, INTENT(OUT)   :: pdf(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins)
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
            ENDIF
          ELSE ! put last point in the last bin
            up = MIN(up,nbins)
            pdf(up) = pdf(up) + C_1_R
          ENDIF
        ENDIF

      ENDDO
    ENDDO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
#endif

  ENDIF

  RETURN
END SUBROUTINE PDF1V2D1G

!########################################################################
! Now the same, but using calculating the PDF over the whole array u
!########################################################################
SUBROUTINE PDF1V3D(ilim, nx,ny,nz, umin_ext,umax_ext,u, nbins,pdf, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx*ny*nz)
  TREAL, INTENT(OUT)   :: pdf(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins)

  ! -------------------------------------------------------------------
  TINTEGER i, up
  TREAL umin, umax, ustep
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
        ENDIF
      ELSE ! put last point in the last bin
        up = MIN(up,nbins)
        pdf(up) = pdf(up) + C_1_R
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

!########################################################################
! And the conditioned one
!########################################################################
SUBROUTINE PDF1V3D1G(ilim, nx,ny,nz, igate,gate, umin_ext,umax_ext,u, nbins,pdf, wrk1d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, nx,ny,nz, nbins
  TREAL umin_ext, umax_ext
  TREAL, INTENT(IN)    :: u(nx*ny*nz)
  TREAL, INTENT(OUT)   :: pdf(nbins+2) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk1d(nbins)
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
          ENDIF
        ELSE ! put last point in the last bin
          up = MIN(up,nbins)
          pdf(up) = pdf(up) + C_1_R
        ENDIF
      ENDIF

    ENDDO

#ifdef USE_MPI
    impi = nbins
    CALL MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins) = wrk1d(1:nbins)
#endif

  ENDIF

  RETURN
END SUBROUTINE PDF1V3D1G

!########################################################################
! Joint PDFs
!########################################################################
SUBROUTINE PDF2V2D(nx,ny,nz, j, u,v, nbins,pdf, wrk2d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, j, nbins(2)
  TREAL, INTENT(IN)    :: u(nx,ny,nz), v(nx,ny,nz)
  TREAL, INTENT(OUT)   :: pdf(nbins(1)*nbins(2) +2 +2*nbins(1)) ! Space at the end for min/max values of sample variable
  TREAL, INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))              ! nbins(2) should be greater than 2 for enough memory space

  ! -------------------------------------------------------------------
  TINTEGER i,k, up,vp, ip, offset
  TREAL umin,umax,ustep

#ifdef USE_MPI
  INTEGER ims_err, impi
  TREAL umin_p, umax_p
#endif

  ! ###################################################################
  pdf = C_0_R

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
  IF ( ustep .EQ. C_0_R) THEN ! Just 1 point, force all in first bin; instead of factor 2, any factor >1 would work
    ustep = (umax-umin) *C_2_R
  ENDIF

  pdf(nbins(1)*nbins(2)+1) = umin + ustep/C_2_R ! Calculate coordinate of histogram
  pdf(nbins(1)*nbins(2)+2) = umax - ustep/C_2_R

  ! Second variable
  DO k = 1,nz
    DO i = 1,nx
      up = INT((u(i,j,k)-umin)/ustep) + 1
      up = MAX(1,MIN(up,nbins(1)))
      ip = offset +up;   pdf(ip) = MIN(pdf(ip),v(i,j,k))
      ip = ip +nbins(1); pdf(ip) = MAX(pdf(ip),v(i,j,k))
    ENDDO
  ENDDO
#ifdef USE_MPI
  impi = nbins(1)
  ip = offset +1;    CALL MPI_ALLREDUCE(pdf(ip), wrk2d, impi, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
  pdf(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
  ip = ip +nbins(1); CALL MPI_ALLREDUCE(pdf(ip), wrk2d, impi, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  pdf(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
#endif

  DO up = 1,nbins(1) ! Calculate Step in Histogram
    wrk2d(up) = ( pdf(offset+up+nbins(1)) -pdf(offset+up) ) /M_REAL(nbins(2))
    IF ( wrk2d(up) .EQ. C_0_R ) THEN ! Just 1 point, force all in first bin
      wrk2d(up) = ( pdf(offset+up+nbins(1)) -pdf(offset+up) ) *C_2_R
    ENDIF
  ENDDO

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
    ENDDO
  ENDDO

  DO up = 1,nbins(1) ! Calculate coordinate in the histogram; I need pdf(offset+up) before
    ip = offset +up;   pdf(ip) = pdf(ip) + wrk2d(up) /C_2_R
    ip = ip +nbins(1); pdf(ip) = pdf(ip) - wrk2d(up) /C_2_R
  ENDDO

#ifdef USE_MPI
  impi = nbins(1)*nbins(2)
  CALL MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
#endif

  RETURN
END SUBROUTINE PDF2V2D

!########################################################################
! Now the same, but using calculating the PDF over the whole array u
!########################################################################
SUBROUTINE PDF2V3D(nx,ny,nz, u,v, nbins,pdf, wrk2d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, nbins(2)
  TREAL, INTENT(IN)    :: u(nx*ny*nz), v(nx*ny*nz)
  TREAL, INTENT(OUT)   :: pdf(nbins(1)*nbins(2) +2 +2*nbins(1)) ! Space at the end for min/max values of sample variable
  TREAL, INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))              ! nbins(2) should be greater than 2 for enough memory space

  ! -------------------------------------------------------------------
  TINTEGER i,k, up,vp, ip, offset
  TREAL umin,umax,ustep
#ifdef USE_MPI
  INTEGER ims_err, impi
#endif

  ! ###################################################################
  pdf = C_0_R

  offset = nbins(1)*nbins(2) +2

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  ! First variable
  CALL MINMAX(nx,ny,nz, u, umin,umax)

  ustep = (umax-umin) /M_REAL(nbins(1)) ! Calculate Step in Histogram
  IF ( ustep .EQ. C_0_R) THEN ! Just 1 point, force all in first bin; instead of factor 2, any factor >1 would work
    ustep = (umax-umin) *C_2_R
  ENDIF

  pdf(nbins(1)*nbins(2)+1) = umin       ! We use umin/umax for code readibility
  pdf(nbins(1)*nbins(2)+2) = umax

  ! Second variable
  DO i = 1,nx*ny*nz
    up = INT((u(i)-umin)/ustep) + 1
    up = MAX(1,MIN(up,nbins(1)))
    ip = offset +up;   pdf(ip) = MIN(pdf(ip),v(i))
    ip = ip +nbins(1); pdf(ip) = MAX(pdf(ip),v(i))
  ENDDO
#ifdef USE_MPI
  impi = nbins(1)
  ip = offset +1;    CALL MPI_ALLREDUCE(pdf(ip), wrk2d, impi, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
  pdf(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
  ip = ip +nbins(1); CALL MPI_ALLREDUCE(pdf(ip), wrk2d, impi, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  pdf(ip:ip+nbins(1)) = wrk2d(1:nbins(1))
#endif

  DO up = 1,nbins(1) ! Calculate Step in Histogram
    wrk2d(up) = ( pdf(offset+up+nbins(1)) -pdf(offset+up) ) /M_REAL(nbins(2))
    IF ( wrk2d(up) .EQ. C_0_R ) THEN ! Just 1 point, force all in first bin
      wrk2d(up) = ( pdf(offset+up+nbins(1)) -pdf(offset+up) ) *C_2_R
    ENDIF
  ENDDO

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  DO i = 1,nx*ny*nz
    up = INT((u(i)-umin)          /ustep    ) + 1
    up = MAX(1,MIN(up,nbins(1)))
    vp = INT((v(i)-pdf(offset+up))/wrk2d(up)) + 1
    vp = MAX(1,MIN(vp,nbins(2)))
    ip = (vp-1)*nbins(1) +up
    pdf(ip) = pdf(ip) + C_1_R
  ENDDO

  DO up = 1,nbins(1) ! Calculate coordinate in the histogram; I need pdf(offset+up) before
    ip = offset +up;   pdf(ip) = pdf(ip) + wrk2d(up) /C_2_R
    ip = ip +nbins(1); pdf(ip) = pdf(ip) - wrk2d(up) /C_2_R
  ENDDO

#ifdef USE_MPI
  impi = nbins(1)*nbins(2)
  CALL MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
#endif

  RETURN
END SUBROUTINE PDF2V3D
