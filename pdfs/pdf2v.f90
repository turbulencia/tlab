#include "types.h"

SUBROUTINE PDF2V2D(nx,ny,nz, j, u,v, nbins, pdf, wrk2d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, j, nbins(2)
  TREAL, INTENT(IN)    :: u(nx,ny,nz), v(nx,ny,nz)
  TREAL, INTENT(OUT)   :: pdf(nbins(1)*nbins(2)+4) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))

  ! -------------------------------------------------------------------
  TINTEGER i,k, up,vp, ip
  TREAL umin,umax,ustep,  vmin,vmax,vstep
  
#ifdef USE_MPI
  INTEGER ims_err, impi
  TREAL umin_p, umax_p, vmin_p, vmax_p
#endif

  ! ###################################################################
  pdf = C_0_R

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  umin = u(1,j,1)
  umax = u(1,j,1)
  vmin = v(1,j,1)
  vmax = v(1,j,1)
  DO k = 1,nz
    DO i = 1,nx
      umin = MIN(umin, u(i,j,k))
      umax = MAX(umax, u(i,j,k))
      vmin = MIN(vmin, v(i,j,k))
      vmax = MAX(vmax, v(i,j,k))
    ENDDO
  ENDDO

#ifdef USE_MPI
  CALL MPI_ALLREDUCE(vmin, vmin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
  vmin = vmin_p
  CALL MPI_ALLREDUCE(vmax, vmax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  vmax = vmax_p
  CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
  umin = umin_p
  CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  umax = umax_p
#endif

  ! Calculate Step in Histogram
  ustep = (umax-umin)/M_REAL(nbins(1))
  vstep = (vmax-vmin)/M_REAL(nbins(2))

  ! Calculate x,y coordinate of histogram
  pdf(nbins(1)*nbins(2)+1) = umin + ustep/C_2_R
  pdf(nbins(1)*nbins(2)+2) = umax - ustep/C_2_R
  pdf(nbins(1)*nbins(2)+3) = vmin + vstep/C_2_R
  pdf(nbins(1)*nbins(2)+4) = vmax - vstep/C_2_R

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  IF ( ABS(ustep) .GT. C_0_R .AND. ABS(vstep) .GT. C_0_R ) THEN
    DO k = 1,nz
      DO i = 1,nx
        up = INT((u(i,j,k)-umin)/ustep) + 1
        up = MAX(1,MIN(up,nbins(1)))
        vp = INT((v(i,j,k)-vmin)/vstep) + 1
        vp = MAX(1,MIN(vp,nbins(2)))
        ip = (vp-1)*nbins(1) +up
        pdf(ip) = pdf(ip) + C_1_R
      ENDDO
    ENDDO

#ifdef USE_MPI
    impi = nbins(1)*nbins(2)
    CALL MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
#endif

  ENDIF

  RETURN
END SUBROUTINE PDF2V2D

!########################################################################
! Now the same, but using calculating the PDF over the whole array u
!########################################################################
SUBROUTINE PDF2V3D(nx,ny,nz, u,v, nbins, pdf, wrk2d)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, nbins(2)
  TREAL, INTENT(IN)    :: u(nx*ny*nz), v(nx*ny*nz)
  TREAL, INTENT(OUT)   :: pdf(nbins(1)*nbins(2)+4) ! Space at the end for the min and max values in the sample variable
  TREAL, INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))

  ! -------------------------------------------------------------------
  TINTEGER i, up,vp, ip
  TREAL umin,umax,ustep, vmin,vmax,vstep
#ifdef USE_MPI
  INTEGER ims_err, impi
#endif

  ! ###################################################################
  pdf = C_0_R

  ! -------------------------------------------------------------------
  ! Calculate Minimum and Maximum
  ! -------------------------------------------------------------------
  CALL MINMAX(nx,ny,nz, u, umin,umax)
  CALL MINMAX(nx,ny,nz, v, vmin,vmax)

  ! Calculate Step in Histogram
  ustep = (umax-umin)/M_REAL(nbins(1))
  vstep = (vmax-vmin)/M_REAL(nbins(2))

  ! Calculate x,y coordinate of histogram
  pdf(nbins(1)*nbins(2)+1) = umin + ustep/C_2_R
  pdf(nbins(1)*nbins(2)+2) = umax - ustep/C_2_R
  pdf(nbins(1)*nbins(2)+3) = vmin + vstep/C_2_R
  pdf(nbins(1)*nbins(2)+4) = vmax - vstep/C_2_R

  ! -------------------------------------------------------------------
  ! Calculate Histogram
  ! -------------------------------------------------------------------
  IF ( ABS(ustep) .GT. C_0_R .AND. ABS(vstep) .GT. C_0_R ) THEN
    DO i = 1,nx*ny*nz
      up = INT((u(i)-umin)/ustep) + 1
      up = MAX(1,MIN(up,nbins(1)))
      vp = INT((v(i)-vmin)/vstep) + 1
      vp = MAX(1,MIN(vp,nbins(2)))
      ip = (vp-1)*nbins(1) +up
      pdf(ip) = pdf(ip) + C_1_R
    ENDDO

#ifdef USE_MPI
    impi = nbins(1)*nbins(2)
    CALL MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2))
#endif

  ENDIF

  RETURN
END SUBROUTINE PDF2V3D
