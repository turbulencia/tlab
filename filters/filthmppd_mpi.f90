      SUBROUTINE FILTHMPPD_MPI(kmax, ijmax, nx,  z1, zf1, wrk)

! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter
! # Midpoint(trapezoidal) rule
! # Periodic boundary conditions
! # PEs explicit version to avoid transpose by MPI
! # nx must be less or equal than kmax !
! #
! #########################################################

      IMPLICIT NONE

#include "types.h"
#include "integers.h"

      TINTEGER kmax, ijmax
      TINTEGER nx
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)
      TREAL wrk(ijmax,*)

      TREAL dmul05, dmul1
      TINTEGER k,ij
      TINTEGER ii

      dmul1  = C_1_R/M_REAL(nx)
      dmul05 = C_05_R/M_REAL(nx)

! filling array
      DO k = 1,kmax
         DO ij = 1,ijmax
            wrk(ij,k+nx/2) = z1(ij,k)
         ENDDO
      ENDDO

      DO k = 1,kmax
! first point
         ii = k
         DO ij = 1,ijmax
            zf1(ij,k) = dmul05*wrk(ij,ii)
         ENDDO
! inner points
         DO ii = k+1,k+nx-1
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+dmul1*wrk(ij,ii)
            ENDDO
         ENDDO
! last point
         ii = k+nx
         DO ij = 1,ijmax
            zf1(ij,k) = zf1(ij,k)+dmul05*wrk(ij,ii)
         ENDDO
      ENDDO

      RETURN
      END
