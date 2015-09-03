#include "types.h"

!########################################################################
! Adding in k of the matrix a for all the elements in j
!########################################################################
SUBROUTINE SUM1V1D_V(ny,nz, a, avg,wrk)
  
  IMPLICIT NONE
  
#ifdef USE_MPI
#include "mpif.h"
#endif
  
  TINTEGER,                INTENT(IN)    :: ny,nz
  TREAL, DIMENSION(ny,nz), INTENT(IN)    :: a
  TREAL, DIMENSION(ny),    INTENT(OUT)   :: avg
  TREAL, DIMENSION(ny),    INTENT(INOUT) :: wrk

! -------------------------------------------------------------------
#ifdef USE_MPI
  INTEGER ims_err, len
#endif

  TINTEGER j, k

! ###################################################################
  avg = C_0_R

  DO k = 1,nz
     DO j = 1,ny
        avg(j) = avg(j) + a(j,k)
     ENDDO
  ENDDO

#ifdef USE_MPI
  len = ny
  CALL MPI_ALLREDUCE(avg, wrk, len, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  avg = wrk
#endif

  RETURN
END SUBROUTINE SUM1V1D_V
