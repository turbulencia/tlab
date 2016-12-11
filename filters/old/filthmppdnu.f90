SUBROUTINE FILTHMPPDNU(kmax, ijmax, nx, cf, z1, zf1)

  ! #########################################################
  ! # FILTER LIBRARY
  ! #
  ! # Top-hat filter on a nonuniform mesh
  ! # Midpoint(trapezoidal) rule
  ! # Periodic boundary conditions
  ! # 
  ! #########################################################

  IMPLICIT NONE

#include "types.h"
#include "integers.h"

  TINTEGER kmax, ijmax
  TINTEGER nx
  TREAL cf(nx+1,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER k,ij,ic
  TINTEGER ii,im

  ! ##################################################
  ! # First nx/2 points where periodicity is imposed #
  ! ##################################################

  DO k = 1,nx/2
     ii = k-nx/2
     im = ii+kmax-i1
     im = MOD(im,kmax)+i1
     ic = ii-k+nx/2+1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,im)*cf(ic,k)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        im = ii+kmax-i1
        im = MOD(im,kmax)+i1
        ic = ii-k+nx/2+1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,im)*cf(ic,k)
        ENDDO
     ENDDO
  ENDDO

  ! #################################################
  ! # Interior points where boundary does not enter #
  ! #################################################

  DO k = 1+nx/2,kmax-nx/2
     ii = k-nx/2
     ic = ii-k+nx/2+1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*cf(ic,k)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        ic = ii-k+nx/2+1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*cf(ic,k)
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
     ic = ii-k+nx/2+1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,im)*cf(ic,k)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        im = ii+kmax-i1
        im = MOD(im,kmax)+i1
        ic = ii-k+nx/2+1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,im)*cf(ic,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FILTHMPPDNU
