      SUBROUTINE FILTHMPPD(kmax, ijmax, nx,  z1, zf1)

! #########################################################
! # FILTER LIBRARY
! #
! # Top-hat filter
! # Midpoint(trapezoidal) rule
! # Periodic boundary conditions
! # 
! #########################################################

      IMPLICIT NONE

#include "types.h"

      TINTEGER kmax, ijmax
      TINTEGER nx
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)

      TREAL dmul05 
      TREAL dmul1
      TINTEGER i1, i2
      TINTEGER k,ij
      TINTEGER ii,im

      i1 = 1
      i2 = 2
      dmul1  = C_1_R/M_REAL(nx)
      dmul05 = C_05_R/M_REAL(nx)

! ##################################################
! # First nx/2 points where periodicity is imposed #
! ##################################################

      DO k = 1,nx/2
! first point and last point
         ii =  k-nx/2
         ii = ii+kmax-i1
         ii = MOD(ii,kmax)+i1
         im = k+nx/2
         DO ij = 1,ijmax
!        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
            zf1(ij,k) = dmul05*(z1(ij,ii)+z1(ij,im))
         ENDDO
! inner points         
         DO ii = k-nx/2+1,k+nx/2-1
            im = ii+kmax-i1
            im = MOD(im,kmax)+i1
            DO ij = 1,ijmax
!           zf1(ij,k) = dmul05*(zf1(ij,k)+C_2_R*z1(ij,im))
               zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,im)
            ENDDO
         ENDDO
      ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################

      DO k = 1+nx/2,kmax-nx/2
! first point and last point
         ii = k-nx/2
         im = k+nx/2
         DO ij = 1,ijmax
!        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
            zf1(ij,k) = dmul05*(z1(ij,ii)+z1(ij,im))
         ENDDO
! inner points
         DO ii = k-nx/2+1,k+nx/2-1
            DO ij = 1,ijmax
!           zf1(ij,k) = dmul05*(zf1(ij,k)+C_2_R*z1(ij,ii))
               zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
            ENDDO
         ENDDO
      ENDDO
      
! #################################################
! # Last nx/2 points where periodicity is imposed #
! #################################################

      DO k = kmax-nx/2+1,kmax
! first point and last point
         ii = k-nx/2
         im = k+nx/2
         im = im+kmax-i1
         im = MOD(im,kmax)+i1
         DO ij = 1,ijmax
!        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
            zf1(ij,k) = dmul05*(z1(ij,ii)+z1(ij,im))
         ENDDO
! inner points
         DO ii = k-nx/2+1,k+nx/2-1
            im = ii+kmax-i1
            im = MOD(im,kmax)+i1
            DO ij = 1,ijmax
!           zf1(ij,k) = dmul05*(zf1(ij,k)+C_2_R*z1(ij,im))
               zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,im)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
