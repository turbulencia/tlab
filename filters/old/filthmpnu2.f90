SUBROUTINE FILTHMPNU2(kmax, ijmax, cf, z1, zf1)
  
! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter on a nonuniform mesh
! # Midpoint(trapezoidal) rule
! # Free boundary conditions
! # Particularized for a stencil of 3 points, nx=2
! #
! #########################################################
  
  IMPLICIT NONE

#include "types.h"

  TINTEGER kmax, ijmax
  TREAL cf(3,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER ij, k
  TREAL a1, a2, a3

! ######################################
! # First point where bc are imposed #
! ######################################
  a2 = cf(2,1)
  a3 = cf(3,1)
  DO ij = 1,ijmax
     zf1(ij,1) = z1(ij,1)*a2 + z1(ij,2)*a3
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 2,kmax-1
     a1 = cf(1,k)
     a2 = cf(2,k)
     a3 = cf(3,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,k-1)*a1 + z1(ij,k)*a2 + z1(ij,k+1)*a3
     ENDDO
  ENDDO

! #####################################
! # Last point where bc are imposed #
! #####################################
  a1 = cf(1,kmax)
  a2 = cf(2,kmax)
  DO ij = 1,ijmax
     zf1(ij,kmax) = z1(ij,kmax-1)*a1 + z1(ij,kmax)*a2
  ENDDO

  RETURN
END SUBROUTINE FILTHMPNU2
