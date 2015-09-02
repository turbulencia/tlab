SUBROUTINE FILTHMPPD4(kmax, ijmax, z1, zf1)

! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter
! # Midpoint(trapezoidal) rule
! # Periodic boundary conditions
! # Particularized for a stencil of 5 points, nx=4
! # 
! #########################################################
  
  IMPLICIT NONE

#include "types.h"

  TINTEGER kmax, ijmax
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER k,ij

! ###############################################
! # First 2 points where periodicity is imposed #
! ###############################################
  DO ij = 1,ijmax
     zf1(ij,1) = C_0125_R*(z1(ij,kmax-1) + z1(ij,3)) + C_025_R*(z1(ij,kmax) + z1(ij,1) + z1(ij,2))
  ENDDO

  DO ij = 1,ijmax
     zf1(ij,2) = C_0125_R*(z1(ij,kmax) + z1(ij,4)) + C_025_R*(z1(ij,1) + z1(ij,2) + z1(ij,3))
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################

  DO k = 3,kmax-2
     DO ij = 1,ijmax
        zf1(ij,k) = C_0125_R*(z1(ij,k-2) + z1(ij,k+2)) + C_025_R*(z1(ij,k-1)+z1(ij,k)+z1(ij,k+1))
     ENDDO
  ENDDO

! ##############################################
! # Last 2 points where periodicity is imposed #
! ##############################################
  DO ij = 1,ijmax
     zf1(ij,kmax-1) = C_0125_R*(z1(ij,kmax-3) + z1(ij,1)) + &
          C_025_R*(z1(ij,kmax-2) + z1(ij,kmax-1) + z1(ij,kmax))
  ENDDO

  DO ij = 1,ijmax
     zf1(ij,kmax) = C_0125_R*(z1(ij,kmax-2) + z1(ij,2)) + &
          C_025_R*(z1(ij,kmax-1) + z1(ij,kmax) + z1(ij,1))
  ENDDO

  RETURN
END SUBROUTINE FILTHMPPD4
