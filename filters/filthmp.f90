      SUBROUTINE FILTHMP(kmax, ijmax, nx, cf, z1, zf1)

! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter
! # Midpoint(trapezoidal) rule
! # Free boundary conditions
! # 
! #########################################################

      IMPLICIT NONE

#include "types.h"

      TINTEGER kmax, ijmax
      TINTEGER nx
      TREAL cf(2,nx/2)
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)

      TREAL r0, dmul05, dmul1
      TINTEGER i1, i2
      TINTEGER k,ij
      TINTEGER ii,im

      i1 = 1
      i2 = 2
      r0 = C_0_R

      dmul1  = C_1_R/M_REAL(nx)
      dmul05 = C_05_R/M_REAL(nx)

! ##################################################
! # First nx/2 points where periodicity is imposed #
! ##################################################

      DO k = 1,nx/2
! first 2 points reflect bc
         DO ij = 1,ijmax
            zf1(ij,k) = dmul1*(z1(ij,1)*cf(1,k)+z1(ij,2)*cf(2,k))
         ENDDO
! inner points
         DO ii = 3,k+nx/2-1
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
            ENDDO
         ENDDO
! last point
         IF ( nx .GT. 2 ) THEN
            ii =  k+nx/2
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+dmul05*z1(ij,ii)
            ENDDO
         ENDIF
      ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################

      DO k = 1+nx/2,kmax-nx/2
! first point
         ii =  k-nx/2
         DO ij = 1,ijmax
            zf1(ij,k) = dmul05*z1(ij,ii)
         ENDDO
! inner points
         DO ii = k-nx/2+1,k+nx/2-1
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
            ENDDO
         ENDDO
! last point
         ii =  k+nx/2
         DO ij = 1,ijmax
            zf1(ij,k) = zf1(ij,k)+dmul05*z1(ij,ii)
         ENDDO
      ENDDO
      
! #################################################
! # Last nx/2 points where periodicity is imposed #
! #################################################

      DO k = kmax-nx/2+1,kmax
! last 2 points account for the bc
         im = kmax-k+i1
         DO ij = 1,ijmax
            zf1(ij,k) = dmul1*(z1(ij,kmax)*cf(1,im)+&
                 z1(ij,kmax-1)*cf(2,im))
         ENDDO
! inner points
         DO ii = k-nx/2+1,kmax-2
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
            ENDDO
         ENDDO
! first point
         IF ( nx .GT. 2 ) THEN
            ii =  k-nx/2
            DO ij = 1,ijmax
               zf1(ij,k) = zf1(ij,k)+dmul05*z1(ij,ii)
            ENDDO
         ENDIF

      ENDDO

      RETURN
      END
