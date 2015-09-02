SUBROUTINE FILTHMPPD2(kmax, ijmax, z1, zf1)
  
! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter
! # Midpoint(trapezoidal) rule
! # Periodic boundary conditions
! # Particularized for a stencil of 3 points, nx=2
! # 
! #########################################################
  
  IMPLICIT NONE

#include "types.h"

  TINTEGER kmax, ijmax
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER k,ij

! ############################################
! # First point where periodicity is imposed #
! ############################################
  DO ij = 1,ijmax
     zf1(ij,1) = C_025_R*(z1(ij,kmax) + z1(ij,2)) + C_05_R*z1(ij,1)
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 2,kmax-1
     DO ij = 1,ijmax
        zf1(ij,k) = C_025_R*(z1(ij,k-1) + z1(ij,k+1)) + C_05_R*z1(ij,k)
     ENDDO
  ENDDO

! ###########################################
! # Last point where periodicity is imposed #
! ###########################################
  DO ij = 1,ijmax
     zf1(ij,kmax) = C_025_R*(z1(ij,kmax-1) + z1(ij,1)) + C_05_R*z1(ij,kmax)
  ENDDO

  RETURN
END SUBROUTINE FILTHMPPD2
