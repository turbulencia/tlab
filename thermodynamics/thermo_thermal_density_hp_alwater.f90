#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/04/07 - Alberto de Lozar
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate rho from h, p and composition using thermal equation of state.
!# Only valid for the incompressible case (p=cte) and mixture = airwater
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_THERMAL_DENSITY_HP_ALWATER(nx,ny,nz, al_q,al_h,p,rho)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL, ONLY : imixture, WGHT_INV, THERMO_AI

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: al_h
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: al_q
  TREAL,			INTENT(IN)  :: p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho


! -------------------------------------------------------------------
  TINTEGER ij
  TREAL WMEAN_INV
  TREAL T_al !temperature
  TREAL dummy
  

! ###################################################################
! mixture MIXT_TYPE_AIRWATER. al_q(1,2) contains liquid mass fraction
! ###################################################################

  IF ( (imixture .EQ. MIXT_TYPE_AIRWATER) .OR. (imixture .EQ. MIXT_TYPE_SUPSAT) ) THEN !Alberto: everything equally valid for Super Saturation

!$omp parallel default( shared ) private (T_al,WMEAN_INV,ij,dummy)
     dummy = (WGHT_INV(1)-WGHT_INV(2))
!$omp do
     DO ij = 1,nx*ny*nz
        T_al = (al_h(ij) - al_q(ij,2)*THERMO_AI(6,1,3) )/( (1-al_q(ij,1))*THERMO_AI(1,1,2) + (al_q(ij,1)-al_q(ij,2))*THERMO_AI(1,1,1)&
               + al_q(ij,2)* THERMO_AI(1,1,3) )
        WMEAN_INV = WGHT_INV(2) + al_q(ij,1)*dummy - al_q(ij,2)*WGHT_INV(1)
        rho(ij) = p/(WMEAN_INV*T_al)
     ENDDO
!$omp end do
!$omp end parallel

  !WRITE(*,*) 'Temperature' T_al

! ###################################################################
! ONLY INTENDED TP WORK WITH AIRWATER
! ###################################################################
  ELSE 
  
  ENDIF

  RETURN
END SUBROUTINE THERMO_THERMAL_DENSITY_HP_ALWATER
