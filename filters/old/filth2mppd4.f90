      SUBROUTINE FILTH2MPPD4(kmax, ijmax, cf, z1, zf1)

! #########################################################
! # FILTER LIBRARY
! #
! # Double top-hat filter
! # Midpoint(trapezoidal) rule
! # Periodic boundary conditions
! # Particularized for a stencil of 5 points, nx=nx0+nx1=4
! # 
! #########################################################

      IMPLICIT NONE

#include "types.h"

      TINTEGER kmax, ijmax
      TREAL cf(5)
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)

      TINTEGER k, ij
      TREAL a1, a2, a3, a4, a5

      a1 = cf(1)
      a2 = cf(2)
      a3 = cf(3)
      a4 = cf(4)
      a5 = cf(5)

! ###############################################
! # First 2 points where periodicity is imposed #
! ###############################################

      DO ij = 1,ijmax
         zf1(ij,1) = z1(ij,kmax-1)*a1 + z1(ij,kmax)*a2 + z1(ij,1)*a3 +&
              z1(ij,2)*a4 + z1(ij,3)*a5
      ENDDO

      DO ij = 1,ijmax
         zf1(ij,2) = z1(ij,kmax)*a1 + z1(ij,1)*a2 + z1(ij,2)*a3 +&
              z1(ij,3)*a4 + z1(ij,4)*a5
      ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################

      DO k = 3,kmax-2
         DO ij = 1,ijmax
            zf1(ij,k) = z1(ij,k-2)*a1 + z1(ij,k-1)*a2 + z1(ij,k)*a3 +&
                 z1(ij,k+1)*a4 + z1(ij,k+2)*a5
         ENDDO
      ENDDO
      
! ##############################################
! # Last 2 points where periodicity is imposed #
! ##############################################

      DO ij = 1,ijmax
         zf1(ij,kmax-1) = z1(ij,kmax-3)*a1 + z1(ij,kmax-2)*a2 + &
              z1(ij,kmax-1)*a3 + z1(ij,kmax)*a4 + z1(ij,1)*a5
      ENDDO

      DO ij = 1,ijmax
         zf1(ij,kmax) = z1(ij,kmax-2)*a1 + z1(ij,kmax-1)*a2 + &
              z1(ij,kmax)*a3 + z1(ij,1)*a4 + z1(ij,2)*a5
      ENDDO

      RETURN
      END
