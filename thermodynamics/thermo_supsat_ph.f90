#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/11/08 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculating the equilibrium T and for given enthalpy and pressure, and ql
!# It is not iterative since we do not ned the equilibrium state
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_SUPSAT_PH(nx, ny, nz, Temperature, s_in, qsat_vapor, p, h) 

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : THERMO_AI, THERMO_PSAT, NPSAT, WGHT_INV, dsmooth

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz,2), INTENT(IN)  :: s_in !Array with total water and liquid water
  TREAL, INTENT(IN)                         ::  p !Alberto pressure whcis is now is now an scalar  
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: h  !array with the entahlpy
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: qsat_vapor, Temperature  !Exit array with the saturation pressure vapor and the temperature
  
  
  TREAL q_latent, cp_d, cp_v, cl, psat, psat_coeff(20),  rd_ov_rv !Auxialliary variables to spped up
! -------------------------------------------------------------------
  TINTEGER ij, ipsat, dum_npsat
  TREAL t_loc  



!$omp parallel default( none ) private( ij, t_loc, psat, ipsat, cp_v, cp_d, cl, rd_ov_rv, q_latent, psat_coeff,dum_npsat) shared( Temperature, qsat_vapor,s_in, h, nx, ny, nz, THERMO_AI, WGHT_INV, NPSAT, THERMO_PSAT, p)
  

! ###################################################################

  
  q_latent = -THERMO_AI(6,1,3)
  cp_d = THERMO_AI(1,1,2)
  cp_v = THERMO_AI(1,1,1)
  cl = THERMO_AI(1,1,3)
  rd_ov_rv = C_1_R/WGHT_INV(1)
  dum_npsat = NPSAT
  DO ipsat = dum_npsat,1,-1
     psat_coeff(ipsat) = THERMO_PSAT(ipsat)/p  !Assuming pressure defined the code and in the definiction of PSAT in HPa 
  ENDDO
  

!$omp do
  DO ij = 1,nx*ny*nz
     t_loc = (h(ij) + s_in(ij,2)*q_latent )/( (1-s_in(ij,1))*cp_d + (s_in(ij,1)-s_in(ij,2))*cp_v + s_in(ij,2)* cl )  !Temperature                           
     psat = C_0_R                             
     DO ipsat = dum_npsat,1,-1 !Saturation pressure
	psat = psat*t_loc + psat_coeff(ipsat) !psat is scaled witrh p
     ENDDO     
     qsat_vapor(ij) = (C_1_R-s_in(ij,1))*rd_ov_rv/(C_1_R/psat-C_1_R)
     Temperature(ij) = t_loc

   ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE THERMO_SUPSAT_PH
