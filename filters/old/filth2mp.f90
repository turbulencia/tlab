SUBROUTINE FILTH2MP(kmax, ijmax, nx0, nx1, cf, z1, zf1)

  ! #########################################################
  ! # FILTER LIBRARY
  ! #
  ! # Double top-hat filter
  ! # Midpoint(trapezoidal) rule
  ! # Free boundary conditions
  ! # 
  ! # nx is nx0+nx1
  ! #########################################################

  IMPLICIT NONE

#include "types.h"

  TINTEGER kmax, ijmax
  TINTEGER nx0, nx1
  TREAL cf(nx0+nx1+1,nx0/2+nx1/2+1)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TREAL r0
  TINTEGER i1, i2
  TINTEGER k,ij,nx
  TINTEGER ii,im,ic

  i1 = 1
  i2 = 2
  r0 = C_0_R
  nx = nx0 + nx1

  ! ##################################################
  ! # First nx/2 points where periodicity is imposed #
  ! ##################################################

  DO k = 1,nx/2
     ii = 1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*cf(ii,k)
     ENDDO
     !     DO ii = 1,k+nx/2
     DO ii = 2,k+nx/2
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*cf(ii,k)
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
        zf1(ij,k) = z1(ij,ii)*cf(ic,nx/2+1)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        ic = ii - k + nx/2 + 1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*cf(ic,nx/2+1)
        ENDDO
     ENDDO
  ENDDO

  ! #################################################
  ! # Last nx/2 points where periodicity is imposed #
  ! #################################################

  DO k = kmax-nx/2+1,kmax
     ii = k-nx/2
     ic = kmax - ii + 1
     im = kmax - k + 1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*cf(ic,im)
     ENDDO
     !     DO ii = k-nx/2,kmax
     DO ii = k-nx/2+1,kmax
        ic = kmax - ii + 1
        im = kmax - k + 1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*cf(ic,im)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FILTH2MP
