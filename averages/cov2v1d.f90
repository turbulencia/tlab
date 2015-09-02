!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 2007/10/26 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the average of the line ij in array a and b. Instead of
!# calling AVG_K, source code is put in here to avoid a function call.
!# Interprocedural in-lining would have done the same, if available. 
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

FUNCTION COV2V1D(nx, ny, nz, i, j, a, b, dz, length)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx, ny, nz, i, j
  TREAL a(nx, ny, nz)
  TREAL b(nx, ny, nz)
  TREAL dz(nz)
  TREAL length

! -------------------------------------------------------------------
  TINTEGER k
  TREAL COV2V1D

#ifdef USE_MPI
  TREAL sum
  INTEGER ierr
#endif

! ###################################################################
  COV2V1D = C_0_R

  DO k = 1, nz
     COV2V1D = COV2V1D + a(i,j,k)*b(i,j,k)*dz(k)
  ENDDO

#ifdef USE_MPI
  sum = COV2V1D/length
  CALL MPI_ALLREDUCE(sum, COV2V1D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  COV2V1D = COV2V1D/length
#endif

  RETURN
END FUNCTION COV2V1D
