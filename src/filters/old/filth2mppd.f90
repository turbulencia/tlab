SUBROUTINE FILTH2MPPD(kmax, ijmax, nx0, nx1, cf, z1, zf1)

  ! #########################################################
  ! # FILTER LIBRARY
  ! #
  ! # Double top-hat filter
  ! # Midpoint(trapezoidal) rule
  ! # Periodic boundary conditions
  ! # 
  ! # nx is nx0+nx1
  ! #########################################################

  IMPLICIT NONE

#include "types.h"

  TINTEGER kmax, ijmax
  TINTEGER nx0, nx1
  TREAL cf(nx0+nx1+1)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER k,ij,nx
  TINTEGER ii,im,ic

  nx = nx0 + nx1

  ! ##################################################
  ! # First nx/2 points where periodicity is imposed #
  ! ##################################################

  DO k = 1,nx/2
     ii = k-nx/2
     im = ii+kmax-i1
     im = MOD(im,kmax)+i1
     ic = ii - k + nx/2 + 1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,im)*cf(ic)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        im = ii+kmax-i1
        im = MOD(im,kmax)+i1
        ic = ii - k + nx/2 + 1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,im)*cf(ic)
        ENDDO
     ENDDO
  ENDDO

  ! #################################################
  ! # Interior points where boundary does not enter #
  ! #################################################

  DO k = 1+nx/2,kmax-nx/2
     ii = k-nx/2
     ic = ii - k + nx/2 + 1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*cf(ic)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        ic = ii - k + nx/2 + 1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*cf(ic)
        ENDDO
     ENDDO
  ENDDO

  ! #################################################
  ! # Last nx/2 points where periodicity is imposed #
  ! #################################################

  DO k = kmax-nx/2+1,kmax
     ii = k-nx/2
     im = ii+kmax-i1
     im = MOD(im,kmax)+i1
     ic = ii - k + nx/2 + 1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,im)*cf(ic)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        im = ii+kmax-i1
        im = MOD(im,kmax)+i1
        ic = ii - k + nx/2 + 1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,im)*cf(ic)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FILTH2MPPD
