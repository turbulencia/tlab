! #############################################################
! # Adding in k of the matrix a for all the elements in j
! #
! # a      = Array of data to calculate the average value of
! # avg    = output array with the mean values
! #
! # 09/22/00 Juan Pedro Mellado
! #############################################################
#include "types.h"

SUBROUTINE SUM_K_V( jmax, kmax, a, avg, wrk )
  
  IMPLICIT NONE
  
#ifdef USE_MPI
#include "mpif.h"
#endif
  
  TINTEGER jmax, kmax
  TREAL, DIMENSION(jmax, kmax) :: a
  TREAL, DIMENSION(jmax)       :: avg, wrk
#ifdef USE_MPI
  INTEGER ims_err, len
#endif

  TINTEGER j, k

! ############################
! # Calculate the mean value #
! ############################
  avg = C_0_R

  DO k = 1,kmax
     DO j = 1,jmax
        avg(j) = avg(j) + a(j,k)
     ENDDO
  ENDDO

#ifdef USE_MPI
  len = jmax
  CALL MPI_ALLREDUCE(avg, wrk, len, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  avg = wrk
#endif

  RETURN
END SUBROUTINE SUM_K_V
