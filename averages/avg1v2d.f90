#include "types.h"

!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 2008/04/10 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# From AVG_IK and AVG_IK_V, to add imom
!#
!########################################################################
FUNCTION AVG1V2D(nx,ny,nz, j, imom, a)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN) :: nx,ny,nz, j
  TINTEGER,                   INTENT(IN) :: imom ! Moment order
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN) :: a

! -------------------------------------------------------------------
  TINTEGER i, k
  TREAL AVG1V2D
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi, norm_mpi
#endif

! ###################################################################
  AVG1V2D = C_0_R
  DO k = 1,nz
     DO i = 1,nx
        AVG1V2D = AVG1V2D + a(i,j,k)**imom
     ENDDO
  ENDDO

#ifdef USE_MPI
  sum_mpi = AVG1V2D
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V2D,  1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  sum_mpi = M_REAL(nx*nz)
  CALL MPI_ALLREDUCE(sum_mpi, norm_mpi, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  AVG1V2D = AVG1V2D/norm_mpi
#else
  AVG1V2D = AVG1V2D/M_REAL(nx*nz)
#endif

  RETURN
END FUNCTION AVG1V2D

! #######################################################################
! #######################################################################
! Vector form
SUBROUTINE AVG1V2D_V(nx,ny,nz, imom, a, avg, wrk)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TINTEGER,                   INTENT(IN)    :: imom ! Moment order
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN)    :: a
  TREAL, DIMENSION(ny),       INTENT(OUT)   :: avg
  TREAL, DIMENSION(ny),       INTENT(INOUT) :: wrk

! -------------------------------------------------------------------
  TINTEGER i,j,k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi, norm_mpi
#endif

! ###################################################################
  avg(:) = C_0_R

  DO k = 1,nz 
     DO j = 1,ny
        DO i = 1,nx
           avg(j) = avg(j) + a(i,j,k)**imom
        ENDDO
     ENDDO
  ENDDO

#ifdef USE_MPI
  wrk = avg
  CALL MPI_ALLREDUCE(wrk, avg, ny, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  sum_mpi = M_REAL(nx*nz)
  CALL MPI_ALLREDUCE(sum_mpi, norm_mpi, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  avg = avg/norm_mpi
#else
  avg = avg/M_REAL(nx*nz)
#endif

END SUBROUTINE AVG1V2D_V

! ###################################################################
! ###################################################################
FUNCTION AVG2V2D(nx,ny,nz, j, a, b)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN) :: nx,ny,nz, j
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN) :: a, b

! -------------------------------------------------------------------
  TINTEGER i, k
  TREAL AVG2V2D
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi, norm_mpi
#endif

! ###################################################################
  AVG2V2D = C_0_R
  DO k = 1,nz
     DO i = 1,nx
        AVG2V2D = AVG2V2D + a(i,j,k)*b(i,j,k)
     ENDDO
  ENDDO

#ifdef USE_MPI
  sum_mpi = AVG2V2D
  CALL MPI_ALLREDUCE(sum_mpi, AVG2V2D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  sum_mpi = M_REAL(nx*nz)
  CALL MPI_ALLREDUCE(sum_mpi, norm_mpi, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  AVG2V2D = AVG2V2D/norm_mpi
#else
  AVG2V2D = AVG2V2D/M_REAL(nx*nz)
#endif

  RETURN
END FUNCTION AVG2V2D

