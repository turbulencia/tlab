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
!# Calculate rho from T, p and composition using thermal equation of state.
!#
!########################################################################
SUBROUTINE THERMO_THERMAL_DENSITY(nx,ny,nz, s,p,T, rho)

  USE THERMO_VARS, ONLY : imixture, MRATIO
  USE THERMO_VARS, ONLY : NSP, WGHT_INV

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: p,T
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho

! -------------------------------------------------------------------
  TINTEGER ij, is
  TREAL WMEAN_INV

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO ij = 1,nx*ny*nz
        rho(ij) = MRATIO*p(ij)/T(ij)
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
!        rho(ij) = MRATIO*p(ij)*WMEAN/T(ij)
!
!     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER. s(1,2) contains liquid mass fraction
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO ij = 1,nx*ny*nz
        WMEAN_INV = WGHT_INV(2) + s(ij,1)*(WGHT_INV(1)-WGHT_INV(2)) - s(ij,2)*WGHT_INV(1)
        rho(ij) = MRATIO*p(ij)/(WMEAN_INV*T(ij))
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

        rho(ij) = MRATIO*p(ij)/(WMEAN_INV*T(ij))

     ENDDO
  ENDIF

  RETURN
END SUBROUTINE THERMO_THERMAL_DENSITY
