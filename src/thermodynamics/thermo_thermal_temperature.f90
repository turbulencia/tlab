#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/05/07 - J.P. Mellado
!#              Created
!# 2007/06/20 - J.P. Mellado
!#              AIRWATER case included.
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate T from rho, p and composition using thermal equation of state.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_THERMAL_TEMPERATURE(nx,ny,nz, s,p,rho, T)

  USE THERMO_GLOBAL, ONLY : imixture, MRATIO
  USE THERMO_GLOBAL, ONLY : NSP, WGHT_INV

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: p,rho
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: T

! -------------------------------------------------------------------
  TINTEGER ij, is
  TREAL WMEAN_INV

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO ij = 1,nx*ny*nz
        T(ij) = MRATIO*p(ij)/rho(ij)
     ENDDO

! ###################################################################
! Mixture defined by a conserved scalar Z, Y_i=f_i(Z)
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_BS          .OR.&
            imixture .EQ. MIXT_TYPE_BSZELDOVICH .OR.&
            imixture .EQ. MIXT_TYPE_QUASIBS    ) THEN
!     DO ij = 1,nx*ny*nz
!
!#define MACRO_ZINPUT s(ij,inb_scal)
!#include "dns_chem_mass.h"
!
!        T(ij) = MRATIO*p(ij)*WMEAN/rho(ij)
!
!     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER. s(1,2) contains liquid mass fraction
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO ij = 1,nx*ny*nz
        WMEAN_INV = WGHT_INV(2) + s(ij,1)*(WGHT_INV(1)-WGHT_INV(2)) - s(ij,2)*WGHT_INV(1)
        T(ij) = MRATIO*p(ij)/(WMEAN_INV*rho(ij))
     ENDDO

! ###################################################################
! General mixture of NSP species. s contains NSP-1 mass fractions
! ###################################################################
  ELSE 
     DO ij = 1,nx*ny*nz

! Compute mean molecular weight. This algorithm reduces the ops to
! 2(NSP-1) +ops and N-1 xops, one less x ops than before
        WMEAN_INV = C_0_R
        DO is = 1,NSP-1
           WMEAN_INV = WMEAN_INV + s(ij,is)*(WGHT_INV(is)-WGHT_INV(NSP))
        ENDDO
        WMEAN_INV = WMEAN_INV + WGHT_INV(NSP)

        T(ij) = MRATIO*p(ij)/(WMEAN_INV*rho(ij))

     ENDDO

  ENDIF

  RETURN
END SUBROUTINE THERMO_THERMAL_TEMPERATURE
