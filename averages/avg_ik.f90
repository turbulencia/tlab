!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/03/09 - J.P. Mellado
!#              Reducing number of ops.
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the average of the plane j in array a. 
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

FUNCTION AVG_IK(nx,ny,nz, j, a, dx,dz, area)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, j
  TREAL, DIMENSION(*)        :: dx,dz
  TREAL, DIMENSION(nx,ny,nz) :: a
  TREAL area

! -------------------------------------------------------------------
  TINTEGER i, k, idsp, kdsp
  TREAL AVG_IK, sum

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
     ENDDO
     AVG_IK = AVG_IK + sum*dz(kdsp+k)
  ENDDO

#ifdef USE_MPI
  sum = AVG_IK/area
  CALL MPI_ALLREDUCE(sum, AVG_IK, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#else
  AVG_IK = AVG_IK/area
#endif

  RETURN
END FUNCTION AVG_IK
