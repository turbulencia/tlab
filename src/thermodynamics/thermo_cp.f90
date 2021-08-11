!########################################################################
!# Tool DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/09 - J.P. Mellado
!#              Created
!# 2007/06/20 - J.P. Mellado
!#              AIRWATER case included.
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate rho from Cp from composition and gama
!#
!# C_p = \gamma/(\gamma-1)*R
!#
!# and nondimensional
!#
!# C_p = \gamma/(\gamma-1)*(R^0/(C_{p,0}W_0))\sum Y_i/W_i
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

SUBROUTINE THERMO_CP(nx,ny,nz, s,gama, cp)

  USE THERMO_VARS, ONLY : imixture, GRATIO
  USE THERMO_VARS, ONLY : NSP, WGHT_INV

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: gama
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: cp

! -------------------------------------------------------------------
  TINTEGER ij, is
  TREAL WMEAN_INV

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO ij = 1,nx*ny*nz
        cp(ij) = C_1_R
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
!        cp(ij) = gama(ij)*GRATIO/((gama(ij)-C_1_R)*WMEAN)
!
!     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER. s(1,2) contains liquid mass fraction
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO ij = 1,nx*ny*nz
        WMEAN_INV = WGHT_INV(2) + s(ij,1)*(WGHT_INV(1)-WGHT_INV(2)) - s(ij,2)*WGHT_INV(1)
        cp(ij) = gama(ij)*GRATIO*WMEAN_INV/(gama(ij)-C_1_R)
     ENDDO

! ###################################################################
! General mixture of NSP species. s contains NSP-1 mass fractions
! ###################################################################
  ELSE 
     DO ij = 1,nx*ny*nz

! Compute mean molecular weight of mayor species
        WMEAN_INV = C_0_R
        DO is = 1,NSP-1
           WMEAN_INV = WMEAN_INV + s(ij,is)*(WGHT_INV(is)-WGHT_INV(NSP))
        ENDDO
        WMEAN_INV = WMEAN_INV + WGHT_INV(NSP)

        cp(ij) = gama(ij)*GRATIO*WMEAN_INV/(gama(ij)-C_1_R)

     ENDDO

  ENDIF

  RETURN
END SUBROUTINE THERMO_CP
