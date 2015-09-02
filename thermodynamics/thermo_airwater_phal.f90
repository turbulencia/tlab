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
!# Calculating the equilibrium T and q_l for given enthalpy and pressure,
!# iterative based on (rho,e,q_i)
!#
!# New created by Alberto to speed up. Only valid for air water . Assumed that HERMO_AI(6,1,1) = HERMO_AI(6,1,2) =0

!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_PHAL(nx, ny, nz, z1, p, h) 

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : THERMO_AI, THERMO_PSAT, NPSAT, WGHT_INV, dsmooth

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx, ny, nz
  TREAL z1(nx*ny*nz,*), h(*) !Alberto
  TREAL p !Alberto p is now an scalar
  TREAL q_latent, cp_d, cp_v, cl, psat, qsat, dpdtsat, psat_coeff(20), dpsat_dT_coeff(20), rd_ov_rv
! -------------------------------------------------------------------
  TINTEGER ij, iter, ipsat, dum_npsat, niter
  TREAL t_loc, z1_loc(2), h_loc, dummy, der1, der2, func,dsmooth_loc,partial_ql_qt,deltaql



!$omp parallel default( none ) private( ij, z1_loc, h_loc, t_loc, psat, qsat, dummy, deltaql, iter, func, dpdtsat, ipsat, der1, der2, dsmooth_loc, cp_v, cp_d, cl, rd_ov_rv, q_latent, psat_coeff, dum_npsat, niter, dpsat_dT_coeff, partial_ql_qt) shared( z1, h, nx, ny, nz, dsmooth, THERMO_AI, WGHT_INV, NPSAT, THERMO_PSAT, p)
  
 

! ###################################################################
  niter = 4
  
  q_latent = -THERMO_AI(6,1,3)
  cp_d = THERMO_AI(1,1,2)
  cp_v = THERMO_AI(1,1,1)
  cl = THERMO_AI(1,1,3)
  rd_ov_rv = C_1_R/WGHT_INV(1)
  dum_npsat = NPSAT
  DO ipsat = dum_npsat,1,-1
     psat_coeff(ipsat) = THERMO_PSAT(ipsat)/p  !Assuming pressure defined the code and in the definiction of PSAT in HPa 
  ENDDO
  DO ipsat = dum_npsat-1,1,-1
     dpsat_dT_coeff(ipsat) = psat_coeff(ipsat+1)*ipsat
  ENDDO


!$omp do
  DO ij = 1,nx*ny*nz
! -------------------------------------------------------------------
! initialize, q_l=0
! -------------------------------------------------------------------
     z1_loc(1) = z1(ij,1)
     z1_loc(2) = C_0_R
     h_loc = h(ij)
     t_loc = h_loc/(cp_d+z1(ij,1)*(cp_v-cp_d)) !Temperature     
     psat = C_0_R                              !Saturation pressure
     DO ipsat = dum_npsat,1,-1
	psat = psat*t_loc + psat_coeff(ipsat)
     ENDDO
     
     dummy = rd_ov_rv/(C_1_R/psat-C_1_R); 
     qsat = (C_1_R-z1_loc(1))*dummy
     z1_loc(2) = z1_loc(1) - qsat; !New liquid concentration
     
! -------------------------------------------------------------------
! iteration
! -------------------------------------------------------------------
     IF (z1_loc(2) .LT. C_0_R) THEN !Unsaturated mixture
	z1_loc(2) = C_0_R     !Reset liquid to zero
        deltaql = C_0_R       !Needed by the smoothing function
     ELSE
	DO iter = 1,niter	
		der1 = cp_d*(1-z1_loc(1)) + cp_v*(z1_loc(1)-z1_loc(2)) + z1_loc(2)*cl !Partial derivative respective to the temperature
		func = der1*t_loc-z1_loc(2)*q_latent-h_loc;  ! Function to minimize
		dpdtsat = C_0_R
		DO ipsat = dum_npsat-1,1,-1
			dpdtsat = dpdtsat*t_loc + dpsat_dT_coeff(ipsat)
		ENDDO        
		der2 = -((cl-cp_v)*t_loc-q_latent)*(C_1_R-z1_loc(1))*rd_ov_rv/((C_1_R-psat)*(C_1_R-psat))*dpdtsat; !Partial derivative respective to ql
		t_loc = t_loc - func/(der1+der2);  ! New guess for  Tnew using Newton (derivative = der1+der2)
		psat = C_0_R          !Recalculate saturation pressure at the new temperature
		DO ipsat = dum_npsat,1,-1
	 		psat = psat*t_loc + psat_coeff(ipsat)
		ENDDO
		dummy = rd_ov_rv/(C_1_R/psat-C_1_R); 
                qsat = (C_1_R-z1_loc(1))*dummy
     		z1_loc(2) = z1_loc(1) - qsat; !New Liquid concentration		
	ENDDO
	 deltaql = -dummy*z1_loc(2) !Needed by the smoothing function
     ENDIF     

     IF ( dsmooth .GT. C_0_R ) THEN
           dsmooth_loc = dsmooth*qsat
	   partial_ql_qt = dummy+C_1_R
           z1(ij,2) = partial_ql_qt*dsmooth_loc*log(exp((z1_loc(1)-qsat)/dsmooth_loc)+C_1_R) + deltaql;
     ELSE
	   z1(ij,2) = z1_loc(2)
     ENDIF

     

   ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE THERMO_AIRWATER_PHAL
