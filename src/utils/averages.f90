#include "types.h"

!########################################################################
! Average along line ij
!########################################################################
TREAL FUNCTION AVG1V1D(nx,ny,nz, i,j, imom, a)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_comm_z, ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN) :: nx,ny,nz, i,j, imom ! Moment order
  TREAL,    INTENT(IN) :: a(nx,ny,nz)

  ! -------------------------------------------------------------------
  TINTEGER k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

  ! ###################################################################
  AVG1V1D = C_0_R
  DO k = 1,nz
    AVG1V1D = AVG1V1D + a(i,j,k)**imom
  END DO

  AVG1V1D = AVG1V1D /M_REAL(nz)
#ifdef USE_MPI
  sum_mpi = AVG1V1D /M_REAL(ims_npro_k)
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V1D,  1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
#endif

  RETURN
END FUNCTION AVG1V1D

!########################################################################
! Adding in k of the matrix a for all the elements in j
!########################################################################
SUBROUTINE SUM1V1D_V(ny,nz, a, avg, wrk)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN   ) :: ny,nz
  TREAL,    INTENT(IN   ) :: a(ny,nz)
  TREAL,    INTENT(  OUT) :: avg(ny)
  TREAL,    INTENT(INOUT) :: wrk(ny)

  ! -------------------------------------------------------------------
  TINTEGER j, k
#ifdef USE_MPI
  INTEGER ims_err
#endif

  ! ###################################################################
  avg = C_0_R
  DO k = 1,nz
    DO j = 1,ny
      avg(j) = avg(j) + a(j,k)
    END DO
  END DO

#ifdef USE_MPI
  wrk = avg
  CALL MPI_ALLREDUCE(wrk, avg, ny, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END SUBROUTINE SUM1V1D_V

!########################################################################
! Covariance along line ij
!########################################################################
TREAL FUNCTION COV2V1D(nx,ny,nz, i,j, a,b)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_comm_z, ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN) :: nx,ny,nz, i,j
  TREAL,    INTENT(IN) :: a(nx,ny,nz), b(nx,ny,nz)

  ! -------------------------------------------------------------------
  TINTEGER k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

  ! ###################################################################
  COV2V1D = C_0_R
  DO k = 1,nz
    COV2V1D = COV2V1D + a(i,j,k) *b(i,j,k)
  END DO

  COV2V1D = COV2V1D /M_REAL(nz)
#ifdef USE_MPI
  sum_mpi = COV2V1D /M_REAL(ims_npro_k)
  CALL MPI_ALLREDUCE(sum_mpi, COV2V1D,  1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
#endif

  RETURN
END FUNCTION COV2V1D

!########################################################################
! Average within the plane j
!########################################################################
TREAL FUNCTION AVG1V2D(nx,ny,nz, j, imom, a)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i,ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN) :: nx,ny,nz, j, imom ! Moment order
  TREAL,    INTENT(IN) :: a(nx,ny,nz)

  ! -------------------------------------------------------------------
  TINTEGER i,k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

  ! ###################################################################
  AVG1V2D = C_0_R
  DO k = 1,nz
    DO i = 1,nx
      AVG1V2D = AVG1V2D + a(i,j,k)**imom
    END DO
  END DO

  AVG1V2D = AVG1V2D /M_REAL(nx*nz)
#ifdef USE_MPI
  sum_mpi = AVG1V2D /M_REAL(ims_npro_i*ims_npro_k)
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V2D,  1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END FUNCTION AVG1V2D

! #######################################################################
! #######################################################################
! Vector form
SUBROUTINE AVG1V2D_V(nx,ny,nz, imom, a, avg, wrk)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i,ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN   ) :: nx,ny,nz, imom ! Moment order
  TREAL,    INTENT(IN   ) :: a(nx,ny,nz)
  TREAL,    INTENT(  OUT) :: avg(ny)
  TREAL,    INTENT(INOUT) :: wrk(ny)

  ! -------------------------------------------------------------------
  TINTEGER i,j,k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

  ! ###################################################################
  avg = C_0_R
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1,nx
        avg(j) = avg(j) + a(i,j,k)**imom
      END DO
    END DO
  END DO

  avg = avg /M_REAL(nx*nz)
#ifdef USE_MPI
  wrk = avg /M_REAL(ims_npro_i*ims_npro_k)
  CALL MPI_ALLREDUCE(wrk, avg, ny, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

END SUBROUTINE AVG1V2D_V

!########################################################################
! Average within the plane j conditioned on the intermittency signal given by array gate
!########################################################################
TREAL FUNCTION AVG1V2D1G(nx,ny,nz, j, igate, imom, a, gate)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,   INTENT(IN) :: nx,ny,nz, j, imom   ! Moment order
  INTEGER(1), INTENT(IN) :: igate               ! Gate level to use
  TREAL,      INTENT(IN) :: a(nx,ny,nz)
  INTEGER(1), INTENT(IN) :: gate(nx,ny,nz)

  ! -------------------------------------------------------------------
  TINTEGER i,k, nsample

#ifdef USE_MPI
  TINTEGER nsample_mpi
  TREAL sum_mpi
  INTEGER ims_err
#endif

  ! ###################################################################
  AVG1V2D1G = C_0_R
  nsample = 0
  DO k = 1,nz
    DO i = 1,nx
      IF ( gate(i,j,k) == igate ) THEN
        AVG1V2D1G = AVG1V2D1G + a(i,j,k)**imom
        nsample = nsample+1
      END IF
    END DO
  END DO

#ifdef USE_MPI
  nsample_mpi = nsample
  CALL MPI_ALLREDUCE(nsample_mpi, nsample, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)

  sum_mpi = AVG1V2D1G
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V2D1G, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  IF ( nsample > 0 ) AVG1V2D1G = AVG1V2D1G /M_REAL(nsample)

  RETURN
END FUNCTION AVG1V2D1G

!########################################################################
! Intermittency factor within the plane j
!########################################################################
TREAL FUNCTION INTER1V2D(nx,ny,nz, j, igate, gate)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i,ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,   INTENT(IN) :: nx,ny,nz, j
  INTEGER(1), INTENT(IN) :: gate(nx,ny,nz), igate

  ! -------------------------------------------------------------------
  TINTEGER i,k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

  ! ###################################################################
  INTER1V2D = C_0_R
  DO k = 1,nz
    DO i = 1,nx
      IF ( gate(i,j,k) == igate ) THEN
        INTER1V2D = INTER1V2D + C_1_R
      END IF
    END DO
  END DO

  INTER1V2D = INTER1V2D /M_REAL(nx*nz)
#ifdef USE_MPI
  sum_mpi = INTER1V2D /M_REAL(ims_npro_i*ims_npro_k)
  CALL MPI_ALLREDUCE(sum_mpi, INTER1V2D,  1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END FUNCTION INTER1V2D

! ###################################################################
! Covariance within the plane j
! ###################################################################
TREAL FUNCTION COV2V2D(nx,ny,nz, j, a, b)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i,ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN) :: nx,ny,nz, j
  TREAL,    INTENT(IN) :: a(nx,ny,nz), b(nx,ny,nz)

  ! -------------------------------------------------------------------
  TINTEGER i,k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

  ! ###################################################################
  COV2V2D = C_0_R
  DO k = 1,nz
    DO i = 1,nx
      COV2V2D = COV2V2D + a(i,j,k)*b(i,j,k)
    END DO
  END DO

  COV2V2D = COV2V2D /M_REAL(nx*nz)
#ifdef USE_MPI
  sum_mpi = COV2V2D /M_REAL(ims_npro_i*ims_npro_k)
  CALL MPI_ALLREDUCE(sum_mpi, COV2V2D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END FUNCTION COV2V2D

!########################################################################
!# DESCRIPTION
!#
!# Calculate the average of the plane j in array a over nonuniform grids.
!#
!########################################################################
TREAL FUNCTION AVG_IK(nx,ny,nz, j, a, dx,dz, area)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN) :: nx,ny,nz, j
  TREAL,    INTENT(IN) :: dx(*),dz(*), area
  TREAL,    INTENT(IN) :: a(nx,ny,nz)

  ! -------------------------------------------------------------------
  TINTEGER i,k, idsp,kdsp
  TREAL sum

  ! ###################################################################
#ifdef USE_MPI
  idsp = ims_offset_i
  kdsp = ims_offset_k
#else
  idsp = 0
  kdsp = 0
#endif

  ! number of + ops: nx*nz*1 + nz*1
  ! number of * ops: nx*nz*1 + nz*1
  AVG_IK = C_0_R
  DO k = 1,nz
    sum = C_0_R
    DO i = 1,nx
      sum = sum + a(i,j,k)*dx(idsp+i)
    END DO
    AVG_IK = AVG_IK + sum*dz(kdsp+k)
  END DO

  AVG_IK = AVG_IK /area
#ifdef USE_MPI
  sum = AVG_IK
  CALL MPI_ALLREDUCE(sum, AVG_IK, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#else
#endif

  RETURN
END FUNCTION AVG_IK

!########################################################################
!########################################################################
! Vector form
SUBROUTINE AVG_IK_V(nx,ny,nz, jm, a, dx,dz, avg, wrk, area)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER, INTENT(IN   ) :: nx,ny,nz, jm
  TREAL,    INTENT(IN   ) :: dx(*),dz(*), area
  TREAL,    INTENT(IN   ) :: a(nx,ny,nz)
  TREAL,    INTENT(  OUT) :: avg(jm)
  TREAL,    INTENT(INOUT) :: wrk(jm)

  ! -------------------------------------------------------------------
  TINTEGER i,j,k, idsp,kdsp
  TREAL sum

  ! ###################################################################
#ifdef USE_MPI
  idsp = ims_offset_i
  kdsp = ims_offset_k
#else
  idsp = 0
  kdsp = 0
#endif

  avg = C_0_R
  DO k = 1,nz
    DO j = 1,jm
      sum = C_0_R
      DO i = 1,nx
        sum = sum + a(i,j,k) *dx(idsp+i)
      END DO
      avg(j) = avg(j) + sum *dz(k+kdsp)
    END DO
  END DO

  avg = avg /area
#ifdef USE_MPI
  wrk = avg
  CALL MPI_ALLREDUCE(wrk, avg, jm, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END SUBROUTINE AVG_IK_V
