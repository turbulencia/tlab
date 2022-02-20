#include "types.h"

!########################################################################
!#
!# Determine Maximum and Minimum Values of a field
!#
!########################################################################
SUBROUTINE MINMAX(imax, jmax, kmax, a, amn, amx)
#ifdef USE_MPI
  USE MPI
#endif

  IMPLICIT NONE

  TINTEGER imax, jmax, kmax
  TREAL a(imax*jmax*kmax), amn, amx

! -----------------------------------------------------------------------
#ifdef USE_MPI
  TREAL pamn, pamx
  INTEGER ims_err
#endif

! #######################################################################
  amn = MINVAL(a)
  amx = MAXVAL(a)

#ifdef USE_MPI
  pamn = amn
  pamx = amx
  CALL MPI_ALLREDUCE(pamn, amn, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
  CALL MPI_ALLREDUCE(pamx, amx, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

  RETURN
END SUBROUTINE MINMAX
