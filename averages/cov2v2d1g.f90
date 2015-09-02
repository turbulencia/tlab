!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 2008/04/10 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the covariance of the plane j between arrays a and b conditioned on the
!# intermittency signal given by array gate
!#
!########################################################################
!# ARGUMENTS 
!#
!# igate    In    Level of the gate signal to use as intermittency function
!# imom     In    Moment order
!#
!########################################################################
FUNCTION COV2V2D1G(nx, ny, nz, j, igate, imom_a, imom_b, a, b, gate, dx, dz)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx, ny, nz, j, imom_a, imom_b
  TREAL dx(nx), dz(nz)
  TREAL a(nx,ny,nz), b(nx,ny,nz)
  INTEGER(1) gate(nx,ny,nz), igate

! -------------------------------------------------------------------
  TINTEGER i, k
  TREAL COV2V2D1G, sum_x, sum_a, area
#ifdef USE_MPI
  INTEGER ims_err
#endif

! ###################################################################
  COV2V2D1G = C_0_R
  area      = C_0_R
  DO k = 1,nz
     sum_a = C_0_R
     sum_x = C_0_R
     DO i = 1,nx
        IF ( gate(i,j,k) .EQ. igate ) THEN 
           sum_a = sum_a + dx(i)*(a(i,j,k)**imom_a)*(b(i,j,k)**imom_b)
           sum_x = sum_x + dx(i)
        ENDIF
     ENDDO
     COV2V2D1G = COV2V2D1G + sum_a*dz(k)
     area      = area      + sum_x*dz(k)
  ENDDO

#ifdef USE_MPI
  sum_x = area
  CALL MPI_ALLREDUCE(sum_x, area, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  sum_a = COV2V2D1G
  CALL MPI_ALLREDUCE(sum_a, COV2V2D1G, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  COV2V2D1G = COV2V2D1G/area

  RETURN
END FUNCTION COV2V2D1G
