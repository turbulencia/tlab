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
!# Calculate thermodynamic properties from h, p and composition in the
!# incompressible formulation, when the thermodynamic pressure is
!# a given profile
!#
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
             ( (C_1_R-q(ij,1))*THERMO_AI(1,1,2) + (q(ij,1)-q(ij,2))*THERMO_AI(1,1,1) &
             + q(ij,2)* THERMO_AI(1,1,3) )
     ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE THERMO_AIRWATER_TEMPERATURE

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_DENSITY(nx,ny,nz, q,h, p, rho)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : WGHT_INV, THERMO_AI, MRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: q   ! total water content, liquid water content
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: h   ! entahlpy
  TREAL, DIMENSION(*),          INTENT(IN)  :: p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho

! -------------------------------------------------------------------
  TINTEGER ij, i, jk
  TREAL WMEAN_INV, P_LOC, T, dummy
  
! ###################################################################
  dummy = (WGHT_INV(1)-WGHT_INV(2))

  ij = 0

  DO jk = 0,ny*nz-1
     P_LOC = MRATIO *p( MOD(jk,ny) +1 )

     DO i = 1,nx
        ij = ij +1

        T = (h(ij) - q(ij,2)*THERMO_AI(6,1,3) ) / &
             ( (C_1_R-q(ij,1))*THERMO_AI(1,1,2) + (q(ij,1)-q(ij,2))*THERMO_AI(1,1,1) &
             + q(ij,2)* THERMO_AI(1,1,3) )
        WMEAN_INV = WGHT_INV(2) + q(ij,1)*dummy - q(ij,2)*WGHT_INV(1)
        rho(ij) = P_LOC /(WMEAN_INV*T)
        
     ENDDO
     
  ENDDO
  
  RETURN
END SUBROUTINE THERMO_AIRWATER_DENSITY

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_BUOYANCY(nx,ny,nz, q,h, p,r, b)

  USE THERMO_GLOBAL, ONLY : WGHT_INV, THERMO_AI, MRATIO
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: q   ! total water content, liquid water content
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: h   ! entahlpy
  TREAL, DIMENSION(*),          INTENT(IN)  :: p,r
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: b

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL WMEAN_INV, P_LOC, R_LOC, R_LOC_INV, T, dummy
  
! ###################################################################
  dummy = (WGHT_INV(1)-WGHT_INV(2))

  ij = 0

  DO jk = 0,ny*nz-1
     is = MOD(jk,ny) +1
     P_LOC = MRATIO *p(is)
     R_LOC = r(is)
     R_LOC_INV = C_1_R /R_LOC
     
     DO i = 1,nx
        ij = ij +1

        T = (h(ij) - q(ij,2)*THERMO_AI(6,1,3) ) / &
             ( (C_1_R-q(ij,1))*THERMO_AI(1,1,2) + (q(ij,1)-q(ij,2))*THERMO_AI(1,1,1) &
             + q(ij,2)* THERMO_AI(1,1,3) )
        WMEAN_INV = WGHT_INV(2) + q(ij,1)*dummy - q(ij,2)*WGHT_INV(1)
        b(ij) = R_LOC_INV *( P_LOC /(WMEAN_INV*T) -R_LOC )
        
     ENDDO
     
  ENDDO
  
  RETURN
END SUBROUTINE THERMO_AIRWATER_BUOYANCY

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_QSAT(nx,ny,nz, q,h, p, T,qsat)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : THERMO_AI, THERMO_PSAT, NPSAT, WGHT_INV

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz,2), INTENT(IN)  :: q       ! total water content, liquid water content
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: h       ! entahlpy
  TREAL, DIMENSION(*),          INTENT(IN)  :: p 
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: qsat, T
  
! -------------------------------------------------------------------
  TINTEGER ij, ipsat, dum_npsat

  TREAL q_latent, cp_d, cp_v, cl, psat, psat_coeff(20),  rd_ov_rv ! Auxialliary variables to spped up
  TREAL t_loc, p_loc

! ###################################################################
  p_loc = p(1) ! To be implemented like in the previous routines
  
!$omp parallel default( none ) &
!$omp private( ij, t_loc, p_loc, psat, ipsat, cp_v, cp_d, cl, rd_ov_rv, q_latent, psat_coeff,dum_npsat) &
!$omp shared( T, qsat,q, h, nx, ny, nz, THERMO_AI, WGHT_INV, NPSAT, THERMO_PSAT, p)

  q_latent = -THERMO_AI(6,1,3)
  cp_d = THERMO_AI(1,1,2)
  cp_v = THERMO_AI(1,1,1)
  cl = THERMO_AI(1,1,3)
  rd_ov_rv = C_1_R/WGHT_INV(1)
  dum_npsat = NPSAT
  DO ipsat = dum_npsat,1,-1
     psat_coeff(ipsat) = THERMO_PSAT(ipsat)/p_loc
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
