SUBROUTINE FILTHMPNU4(kmax, ijmax, cf, z1, zf1)
  
! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter on a nonuniform mesh
! # Midpoint(trapezoidal) rule
! # Free boundary conditions
! # Particularized for a stencil of 5 points, nx=4
! #
! #########################################################
  
  IMPLICIT NONE

#include "types.h"

  TINTEGER kmax, ijmax
  TREAL cf(5,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER ij, k
  TREAL a1, a2, a3, a4, a5

! #########################################
! # First 2 points where bc are imposed #
! #########################################
  a3 = cf(3,1)
  a4 = cf(4,1)
  a5 = cf(5,1)
  DO ij = 1,ijmax
     zf1(ij,1) = z1(ij,1)*a3 + z1(ij,2)*a4 + z1(ij,3)*a5
  ENDDO

  a2 = cf(2,2)
  a3 = cf(3,2)
  a4 = cf(4,2)
  a5 = cf(5,2)
  DO ij = 1,ijmax
     zf1(ij,2) = z1(ij,1)*a2 + z1(ij,2)*a3 + z1(ij,3)*a4 + z1(ij,4)*a5
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 3,kmax-2
     a1 = cf(1,k)
     a2 = cf(2,k)
     a3 = cf(3,k)
     a4 = cf(4,k)
     a5 = cf(5,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,k-2)*a1 + z1(ij,k-1)*a2 + z1(ij,k)*a3 + z1(ij,k+1)*a4 + z1(ij,k+2)*a5
     ENDDO
  ENDDO

! ########################################
! # Last 2 points where bc are imposed #
! ########################################
  a1 = cf(1,kmax-1)
  a2 = cf(2,kmax-1)
  a3 = cf(3,kmax-1)
  a4 = cf(4,kmax-1)
  DO ij = 1,ijmax
     zf1(ij,kmax-1) = z1(ij,kmax-3)*a1 + z1(ij,kmax-2)*a2 + z1(ij,kmax-1)*a3 + z1(ij,kmax)*a4
  ENDDO

  a1 = cf(1,kmax)
  a2 = cf(2,kmax)
  a3 = cf(3,kmax)
  DO ij = 1,ijmax
     zf1(ij,kmax) = z1(ij,kmax-2)*a1 + z1(ij,kmax-1)*a2 + z1(ij,kmax)*a3
  ENDDO

  RETURN
END SUBROUTINE FILTHMPNU4
