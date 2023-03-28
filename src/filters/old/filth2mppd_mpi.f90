      SUBROUTINE FILTH2MPPD_MPI(kmax, ijmax, nx, cf, z1, zf1, wrk)

! #########################################################
! # FILTER LIBRARY
! #
! # Double top-hat filter
! # Midpoint(trapezoidal) rule
! # Periodic boundary conditions
! # PEs Local version to avoid transpose by MPI
! # 
! # nx is nx0+nx1
! #
! #########################################################

      IMPLICIT NONE

#include "types.h"

      TINTEGER kmax, ijmax
      TINTEGER nx
      TREAL cf(nx+1)
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)
      TREAL wrk(ijmax,*)

      TINTEGER k,ij
      TINTEGER ii,ic

! filling array
      DO k = 1,kmax
         DO ij = 1,ijmax
            wrk(ij,k+nx/2) = z1(ij,k)
         ENDDO
      ENDDO

      DO k = 1,kmax
            ii = k
            ic = i1
            DO ij = 1,ijmax
               zf1(ij,k) = wrk(ij,ii)*cf(ic)
            ENDDO
         DO ii = k+1,k+nx
            ic = ii - k + 1
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+wrk(ij,ii)*cf(ic)
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
