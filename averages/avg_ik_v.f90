SUBROUTINE AVG_IK_V(imax, jmax, kmax, jm, a, dx, dz, avg, wrk, area)
!
! This function calculates the average value of the elements
! in the array A on constant J slices
!
! Passed Variables
! ------------------------------------------------------------------------
! imax, jmax, kmax - Dimensions of array
! j - J-station to calculate average of
! a - Array of data to calculate the average value of
! ------------------------------------------------------------------------
!
! *****************************
! *** Variable Declarations ***
! *****************************
!
#include "types.h"

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

!
! *** Passed Variables ***
!
  TINTEGER imax, jmax, kmax, jm
  TREAL a(imax, jmax, kmax)
  TREAL dx(imax), dz(kmax)
  TREAL avg(jm), wrk(jm)
  TREAL area
#ifdef USE_MPI
  INTEGER ims_err, len
#endif

!
! *** Local Variables ***
!
  TINTEGER i, j, k


!
! *******************************
! *** Calculate the RMS value ***
! *******************************
!
  DO j=1,jm
     avg(j) = C_0_R
  ENDDO

  DO k = 1, kmax
     DO j=1, jm
        DO i = 1, imax
           avg(j) = avg(j) + a(i,j,k)*dx(i)*dz(k)
        ENDDO
     ENDDO
  ENDDO

#ifdef USE_MPI
  len = jm
  CALL MPI_REDUCE(avg, wrk, len, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
  DO j=1,jm
     avg(j) = wrk(j)/area
  ENDDO
#else
  DO j=1,jm
     avg(j) = avg(j)/area
  ENDDO
#endif

  RETURN
END SUBROUTINE AVG_IK_V
