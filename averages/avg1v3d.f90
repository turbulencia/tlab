!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 2007/03/09 - J.P. Mellado
!#              Created
!# 2008/04/10 - J.P. Mellado
!#              Moment argument introduced
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the volume average in array a. 
!#
!########################################################################
!# ARGUMENTS 
!#
!# imom     In    Moment order
!#
!########################################################################
#include "types.h"

FUNCTION AVG1V3D(nx,ny,nz, imom, a)

#ifdef USE_MPI
  USE DNS_GLOBAL, ONLY : imax_total, jmax_total, kmax_total
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz, imom
  TREAL a(*)

! -------------------------------------------------------------------
  TINTEGER ij
  TREAL AVG1V3D
#ifdef USE_MPI
  INTEGER ims_err
  TREAL sum_mpi
#endif

! ###################################################################
  AVG1V3D = C_0_R
  DO ij = 1, nx*ny*nz
     AVG1V3D = AVG1V3D + a(ij)**imom
  ENDDO

#ifdef USE_MPI
  sum_mpi = AVG1V3D/M_REAL(imax_total*jmax_total)
  sum_mpi = sum_mpi/M_REAL(kmax_total) ! In two steps is case nx*ny*nz is larger than INT(4)
  CALL MPI_ALLREDUCE(sum_mpi, AVG1V3D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#else
  AVG1V3D = AVG1V3D/M_REAL(nx*ny)
  AVG1V3D = AVG1V3D/M_REAL(nz) ! In two steps is case nx*ny*nz is larger than INT(4)
#endif

  RETURN
END FUNCTION AVG1V3D
