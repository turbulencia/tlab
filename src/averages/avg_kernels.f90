#include "types.h"

!########################################################################
! Average along line ij
!########################################################################
TREAL FUNCTION AVG1V1D(nx,ny,nz, i,j, imom, a)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_comm_z
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, i,j, imom
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN) :: a

! -------------------------------------------------------------------
  TINTEGER k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi, norm_mpi
#endif

! ###################################################################
  AVG1V1D = C_0_R
  DO k = 1, nz
     AVG1V1D = AVG1V1D + a(i,j,k)**imom
  ENDDO

#ifdef USE_MPI
  sum_mpi = AVG1V1D
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V1D,  1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
  sum_mpi = M_REAL(nz)
  CALL MPI_ALLREDUCE(sum_mpi, norm_mpi, 1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
  AVG1V1D = AVG1V1D/norm_mpi
#else
  AVG1V1D = AVG1V1D/M_REAL(nz)
#endif

  RETURN
END FUNCTION AVG1V1D

!########################################################################
! Covariance along line ij
!########################################################################
TREAL FUNCTION COV2V1D(nx,ny,nz, i,j, a,b)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_comm_z
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, i,j
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN) :: a,b

! -------------------------------------------------------------------
  TINTEGER k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi, norm_mpi
#endif

! ###################################################################
  COV2V1D = C_0_R
  DO k = 1, nz
     COV2V1D = COV2V1D + a(i,j,k) *b(i,j,k)
  ENDDO

#ifdef USE_MPI
  sum_mpi = COV2V1D
  CALL MPI_ALLREDUCE(sum_mpi, COV2V1D,  1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
  sum_mpi = M_REAL(nz)
  CALL MPI_ALLREDUCE(sum_mpi, norm_mpi, 1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
  COV2V1D = COV2V1D/norm_mpi
#else
  COV2V1D = COV2V1D/M_REAL(nz)
#endif

  RETURN
END FUNCTION COV2V1D

!########################################################################
! Average within the plane j
!########################################################################
TREAL FUNCTION AVG1V2D(nx,ny,nz, j, imom, a)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_npro_i,ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN) :: nx,ny,nz, j
  TINTEGER,                   INTENT(IN) :: imom ! Moment order
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN) :: a

! -------------------------------------------------------------------
  TINTEGER i,k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi!, norm_mpi
#endif

! ###################################################################
  AVG1V2D = C_0_R
  DO k = 1,nz
     DO i = 1,nx
        AVG1V2D = AVG1V2D + a(i,j,k)**imom
     ENDDO
  ENDDO

  AVG1V2D = AVG1V2D/M_REAL(nx*nz)
#ifdef USE_MPI
  sum_mpi = AVG1V2D/M_REAL(ims_npro_i*ims_npro_k)
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V2D,  1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
! #ifdef USE_MPI
!   sum_mpi = AVG1V2D
!   CALL MPI_ALLREDUCE(sum_mpi, AVG1V2D,  1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
!   sum_mpi = M_REAL(nx*nz)
!   CALL MPI_ALLREDUCE(sum_mpi, norm_mpi, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
!   AVG1V2D = AVG1V2D/norm_mpi
! #else
!   AVG1V2D = AVG1V2D/M_REAL(nx*nz)
! #endif

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

!########################################################################
! Average within the plane j conditioned on the intermittency signal given by array gate
!########################################################################
TREAL FUNCTION AVG1V2D1G(nx,ny,nz, j, igate, imom, a, gate)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                        INTENT(IN) :: nx,ny,nz, j
  INTEGER(1),                      INTENT(IN) :: igate ! Gate level to use
  TINTEGER,                        INTENT(IN) :: imom  ! Moment order
  TREAL, DIMENSION(nx,ny,nz),      INTENT(IN) :: a
  INTEGER(1), DIMENSION(nx,ny,nz), INTENT(IN) :: gate

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
        IF ( gate(i,j,k) .EQ. igate ) THEN
           AVG1V2D1G = AVG1V2D1G + a(i,j,k)**imom
           nsample = nsample+1
        ENDIF
     ENDDO
  ENDDO

#ifdef USE_MPI
  nsample_mpi = nsample
  CALL MPI_ALLREDUCE(nsample_mpi, nsample, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)

  sum_mpi = AVG1V2D1G
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V2D1G, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  IF ( nsample .GT. 0 ) AVG1V2D1G = AVG1V2D1G/M_REAL(nsample)

  RETURN
END FUNCTION AVG1V2D1G

!########################################################################
! Intermittency factor within the plane j 
!########################################################################
TREAL FUNCTION INTER1V2D(nx,ny,nz, j, igate, gate)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_npro_i,ims_npro_k
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
      IF ( gate(i,j,k) .EQ. igate ) THEN
        INTER1V2D = INTER1V2D + C_1_R
      ENDIF
    ENDDO
  ENDDO

  INTER1V2D = INTER1V2D /M_REAL(nx*nz)
#ifdef USE_MPI
  sum_mpi = INTER1V2D/M_REAL(ims_npro_i*ims_npro_k)
  CALL MPI_ALLREDUCE(sum_mpi, INTER1V2D,  1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END FUNCTION INTER1V2D

! ###################################################################
! Covariance within the plane j
! ###################################################################
TREAL FUNCTION COV2V2D(nx,ny,nz, j, a, b)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN) :: nx,ny,nz, j
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN) :: a, b

! -------------------------------------------------------------------
  TINTEGER i,k
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi, norm_mpi
#endif

! ###################################################################
  COV2V2D = C_0_R
  DO k = 1,nz
     DO i = 1,nx
        COV2V2D = COV2V2D + a(i,j,k)*b(i,j,k)
     ENDDO
  ENDDO

#ifdef USE_MPI
  sum_mpi = COV2V2D
  CALL MPI_ALLREDUCE(sum_mpi, COV2V2D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  sum_mpi = M_REAL(nx*nz)
  CALL MPI_ALLREDUCE(sum_mpi, norm_mpi, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  COV2V2D = COV2V2D/norm_mpi
#else
  COV2V2D = COV2V2D/M_REAL(nx*nz)
#endif

  RETURN
END FUNCTION COV2V2D

!########################################################################
! Average of array a
!########################################################################
TREAL FUNCTION AVG1V3D(nx,ny,nz, imom, a)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_npro
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN) :: nx,ny,nz
  TINTEGER,                   INTENT(IN) :: imom ! Moment order
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN) :: a

! -------------------------------------------------------------------
  TINTEGER ij

#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

! ###################################################################
  AVG1V3D = C_0_R
  DO ij = 1, nx*ny*nz
     AVG1V3D = AVG1V3D + a(ij)**imom
  ENDDO

  AVG1V3D = AVG1V3D/M_REAL(nx*ny)
  AVG1V3D = AVG1V3D/M_REAL(nz) ! In two steps is case nx*ny*nz is larger than INT(4)
#ifdef USE_MPI
  sum_mpi = AVG1V3D/M_REAL(ims_npro)
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V3D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END FUNCTION AVG1V3D

!########################################################################
! Average of array a  conditioned on the intermittency signal given by array gate
!########################################################################
TREAL FUNCTION AVG1V3D1G(nx,ny,nz, igate, imom, a, gate)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                        INTENT(IN) :: nx,ny,nz
  INTEGER(1),                      INTENT(IN) :: igate ! Gate level to use
  TINTEGER,                        INTENT(IN) :: imom  ! Moment order
  TREAL, DIMENSION(nx*ny*nz),      INTENT(IN) :: a
  INTEGER(1), DIMENSION(nx*ny*nz), INTENT(IN) :: gate

! -------------------------------------------------------------------
  TINTEGER ij, nsample

#ifdef USE_MPI
  TINTEGER nsample_mpi
  TREAL sum_mpi
  INTEGER ims_err
#endif

! ###################################################################
  AVG1V3D1G = C_0_R
  nsample = 0
  DO ij = 1,nx*ny*nz
     IF ( gate(ij) .EQ. igate ) THEN
        AVG1V3D1G = AVG1V3D1G + a(ij)**imom
        nsample = nsample+1
     ENDIF
  ENDDO

#ifdef USE_MPI
  nsample_mpi = nsample
  CALL MPI_ALLREDUCE(nsample_mpi, nsample, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)

  sum_mpi = AVG1V3D1G
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V3D1G, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  IF ( nsample .GT. 0 ) AVG1V3D1G = AVG1V3D1G/M_REAL(nsample)

  RETURN
END FUNCTION AVG1V3D1G
