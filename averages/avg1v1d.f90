!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/03/09 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the average of the line ij in array a. 
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

FUNCTION AVG1V1D(nx, ny, nz, i, j, a, dz, length)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx, ny, nz, i, j
  TREAL a(nx, ny, nz)
  TREAL dz(nz)
  TREAL length

! -------------------------------------------------------------------
  TINTEGER k
  TREAL AVG1V1D

#ifdef USE_MPI
  TREAL sum
  INTEGER ierr
#endif

! ###################################################################
  AVG1V1D = C_0_R

  DO k = 1, nz
     AVG1V1D = AVG1V1D + a(i,j,k)*dz(k)
  ENDDO

#ifdef USE_MPI
  sum = AVG1V1D/length
  CALL MPI_ALLREDUCE(sum, AVG1V1D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  AVG1V1D = AVG1V1D/length
#endif

  RETURN
END FUNCTION AVG1V1D
