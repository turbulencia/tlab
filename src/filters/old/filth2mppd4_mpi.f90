      SUBROUTINE FILTH2MPPD4_MPI(kmax, ijmax, nx, cf, z1, zf1, wrk)

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
#include "integers.h"

      TINTEGER kmax, ijmax
      TINTEGER nx
      TREAL cf(nx+1)
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)
      TREAL wrk(ijmax,*)

      TINTEGER k,ij

! filling array
      DO k = 1,kmax
         DO ij = 1,ijmax
            wrk(ij,k+nx/2) = z1(ij,k)
         ENDDO
      ENDDO

      DO k = 1,kmax
         DO ij = 1,ijmax
            zf1(ij,k) = wrk(ij,k)*cf(1) + wrk(ij,k+1)*cf(2) +&
                 wrk(ij,k+2)*cf(3) + wrk(ij,k+3)*cf(4) +&
                 wrk(ij,k+4)*cf(5)
         ENDDO
      ENDDO
      
      RETURN
      END
