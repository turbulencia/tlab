#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2007/04/07 - Alberto de Lozar
!#              Created
!# 2016/10/15 - J.P. Mellado
!#              Cleaning
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate rho from h, p and composition for incompressible case.
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_DENSITY(nx,ny,nz, q,p,h, rho)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : WGHT_INV, THERMO_AI

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: q       ! total water content, liquid water content
  TREAL,			INTENT(IN)  :: p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: h       ! entahlpy
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho

! -------------------------------------------------------------------
  TINTEGER ij
  TREAL WMEAN_INV, T, dummy
  
! ###################################################################
! mixture MIXT_TYPE_AIRWATER. q(1,2) contains liquid mass fraction
! ###################################################################
!$omp parallel default( shared ) private (T,WMEAN_INV,ij,dummy)
     dummy = (WGHT_INV(1)-WGHT_INV(2))
!$omp do
     DO ij = 1,nx*ny*nz
        T = (h(ij) - q(ij,2)*THERMO_AI(6,1,3) )/( (1-q(ij,1))*THERMO_AI(1,1,2) + (q(ij,1)-q(ij,2))*THERMO_AI(1,1,1)&
               + q(ij,2)* THERMO_AI(1,1,3) )
        WMEAN_INV = WGHT_INV(2) + q(ij,1)*dummy - q(ij,2)*WGHT_INV(1)
        rho(ij) = p/(WMEAN_INV*T)
     ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE THERMO_AIRWATER_DENSITY
