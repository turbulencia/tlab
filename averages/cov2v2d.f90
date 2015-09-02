!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/11 - J.P. Mellado
!#              Created
!# 2008/04/10 - J.P. Mellado
!#              Moment allowed
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the covariance of the plane j in arrays a and b. Instead of
!# calling AVG_IK, source code is put in here to avoid an array mem call.
!#
!########################################################################
!# ARGUMENTS 
!#
!# imom     In    Moment order
!#
!########################################################################
FUNCTION COV2V2D(nx, ny, nz, j, imom_a, imom_b, a, b, dx, dz, area)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx, ny, nz, j, imom_a, imom_b
  TREAL dx(nx), dz(nz)
  TREAL a(nx, ny, nz), b(nx, ny, nz), area

! -------------------------------------------------------------------
  TINTEGER i, k
  TREAL COV2V2D, sum
#ifdef USE_MPI
  INTEGER ims_err
#endif

! ###################################################################
! number of + ops: nx*nz*1 + nz*1
! number of * ops: nx*nz*2 + nz*1
  COV2V2D = C_0_R
  DO k = 1, nz
     sum = C_0_R
     DO i = 1, nx
        sum = sum + dx(i)*(a(i,j,k)**imom_a)*(b(i,j,k)**imom_b)
     ENDDO
     COV2V2D = COV2V2D + sum*dz(k)
  ENDDO

#ifdef USE_MPI
  sum = COV2V2D/area
  CALL MPI_ALLREDUCE(sum, COV2V2D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#else
  COV2V2D = COV2V2D/area
#endif

  RETURN
END FUNCTION COV2V2D
