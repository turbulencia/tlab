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
!# Calculate thermodynamic properties from h, p and composition.
!#
!# Assumes p is a constant.
!# Assumes that THERMO_AI(6,1,1) = THERMO_AI(6,1,2) = 0
!#
!########################################################################

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_TEMPERATURE(nx,ny,nz, q,h, T)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : THERMO_AI

  IMPLICIT NONE

  TINTEGER, INTENT(IN)                      :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: q   ! total water content, liquid water content
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: h   ! entahlpy
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: T

! -------------------------------------------------------------------
  TINTEGER ij
  
! ###################################################################
!$omp parallel default( shared ) private (ij)
!$omp do
     DO ij = 1,nx*ny*nz
        T(ij) = (h(ij) - q(ij,2)*THERMO_AI(6,1,3) ) / &
             ( (1-q(ij,1))*THERMO_AI(1,1,2) + (q(ij,1)-q(ij,2))*THERMO_AI(1,1,1) &
             + q(ij,2)* THERMO_AI(1,1,3) )
     ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE THERMO_AIRWATER_TEMPERATURE

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_DENSITY(nx,ny,nz, q,p,h, rho)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : WGHT_INV, THERMO_AI

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: q   ! total water content, liquid water content
  TREAL,                        INTENT(IN)  :: p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: h   ! entahlpy
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho

! -------------------------------------------------------------------
  TINTEGER ij
  TREAL WMEAN_INV, T, dummy
  
! ###################################################################
!$omp parallel default( shared ) private (T,WMEAN_INV,ij,dummy)
     dummy = (WGHT_INV(1)-WGHT_INV(2))
!$omp do
     DO ij = 1,nx*ny*nz
        T = (h(ij) - q(ij,2)*THERMO_AI(6,1,3) ) / &
             ( (1-q(ij,1))*THERMO_AI(1,1,2) + (q(ij,1)-q(ij,2))*THERMO_AI(1,1,1) &
             + q(ij,2)* THERMO_AI(1,1,3) )
        WMEAN_INV = WGHT_INV(2) + q(ij,1)*dummy - q(ij,2)*WGHT_INV(1)
        rho(ij) = p/(WMEAN_INV*T)
     ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE THERMO_AIRWATER_DENSITY

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_QSAT(nx,ny,nz, q,p,h, T,qsat)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : THERMO_AI, THERMO_PSAT, NPSAT, WGHT_INV

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz,2), INTENT(IN)  :: q       ! total water content, liquid water content
  TREAL,                        INTENT(IN)  :: p 
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: h       ! entahlpy
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: qsat, T
  
! -------------------------------------------------------------------
  TINTEGER ij, ipsat, dum_npsat

  TREAL q_latent, cp_d, cp_v, cl, psat, psat_coeff(20),  rd_ov_rv ! Auxialliary variables to spped up
  TREAL t_loc  

! ###################################################################
!$omp parallel default( none ) &
!$omp private( ij, t_loc, psat, ipsat, cp_v, cp_d, cl, rd_ov_rv, q_latent, psat_coeff,dum_npsat) &
!$omp shared( T, qsat,q, h, nx, ny, nz, THERMO_AI, WGHT_INV, NPSAT, THERMO_PSAT, p)

  q_latent = -THERMO_AI(6,1,3)
  cp_d = THERMO_AI(1,1,2)
  cp_v = THERMO_AI(1,1,1)
  cl = THERMO_AI(1,1,3)
  rd_ov_rv = C_1_R/WGHT_INV(1)
  dum_npsat = NPSAT
  DO ipsat = dum_npsat,1,-1
     psat_coeff(ipsat) = THERMO_PSAT(ipsat)/p
  ENDDO

!$omp do
  DO ij = 1,nx*ny*nz
     t_loc = (h(ij) + q(ij,2)*q_latent )/( (1-q(ij,1))*cp_d + (q(ij,1)-q(ij,2))*cp_v + q(ij,2)* cl )
     psat = C_0_R                             
     DO ipsat = dum_npsat,1,-1
        psat = psat*t_loc + psat_coeff(ipsat) !psat is scaled with p
     ENDDO     
     qsat(ij) = (C_1_R-q(ij,1))*rd_ov_rv/(C_1_R/psat-C_1_R)
     T(ij) = t_loc

   ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE THERMO_AIRWATER_QSAT
