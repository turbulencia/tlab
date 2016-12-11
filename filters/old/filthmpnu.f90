      SUBROUTINE FILTHMPNU(kmax, ijmax, nx, cf, z1, zf1)

! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter on a nonuniform mesh
! # Midpoint(trapezoidal) rule
! # Free boundary conditions
! # 
! #########################################################

      IMPLICIT NONE

#include "types.h"

      TINTEGER kmax, ijmax
      TINTEGER nx
      TREAL cf(nx+1,kmax)
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)

      TREAL r0, dum
      TINTEGER i1, i2
      TINTEGER k,ij,ic
      TINTEGER ii, iiorg

      i1 = 1
      i2 = 2
      r0 = C_0_R

! ############################################
! # First nx/2 points where bc are imposed #
! ############################################

      DO k = 1,nx/2
         iiorg = nx/2-k+1
            ii = 1
            ic = iiorg + ii
            dum = cf(ic,k)
            DO ij = 1,ijmax
               zf1(ij,k) = z1(ij,ii)*dum
            ENDDO
         DO ii = 2,k+nx/2
            ic = iiorg + ii
            dum = cf(ic,k)
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*dum
            ENDDO
         ENDDO
      ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################

      DO k = 1+nx/2,kmax-nx/2
         iiorg = nx/2-k+1
            ii = k-nx/2
            ic = ii + iiorg
            dum = cf(ic,k)
            DO ij = 1,ijmax
               zf1(ij,k) = z1(ij,ii)*dum
            ENDDO
         DO ii = k-nx/2+1,k+nx/2
            ic = ii + iiorg
            dum = cf(ic,k)
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*dum
            ENDDO
         ENDDO
      ENDDO
      
! ###########################################
! # Last nx/2 points where bc are imposed #
! ###########################################

      DO k = kmax-nx/2+1,kmax
         iiorg = nx/2-k+1
            ii = k-nx/2
            ic = ii + iiorg
            dum = cf(ic,k)
            DO ij = 1,ijmax
               zf1(ij,k) = z1(ij,ii)*dum
            ENDDO
         DO ii = k-nx/2+1,kmax
            ic = ii + iiorg
            dum = cf(ic,k)
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*dum
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
