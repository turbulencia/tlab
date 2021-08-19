#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/21 - J.P. Mellado
!#              Created
!# 2007/07/05 - J.P. Mellado
!#              AirWater case included.
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate gama from T and composition.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_GAMMA(nx,ny,nz, s, T, gama)

  USE THERMO_VARS, ONLY : imixture, gama0, GRATIO
  USE THERMO_VARS, ONLY : NSP, NCP_CHEMKIN, WGHT_INV, THERMO_AI, THERMO_TLIM
  USE THERMO_VARS, ONLY : YMASS

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: T
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: gama

! -------------------------------------------------------------------
  TINTEGER ij, is, icp, im
  TREAL WMEAN_INV, HEAT_CAPACITY, HEAT_CAPACITY_I

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO ij = 1, nx*ny*nz
        gama(ij) = gama0
     ENDDO

! ###################################################################
! Mixture defined by a conserved scalar Z, Y_i=f_i(Z)
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_BS          .OR.&
            imixture .EQ. MIXT_TYPE_BSZELDOVICH .OR.&
            imixture .EQ. MIXT_TYPE_QUASIBS    ) THEN
!     DO ij = 1,nx*ny*nz
!
!#define MACRO_TINPUT T(ij)
!#include "dns_chem_enth.h"
!        gama(ij) = CPW/(CPW-GRATIO)
!
!     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER. s(1,2) contains liquid mass fraction
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO ij = 1, nx*ny*nz
        HEAT_CAPACITY = THERMO_AI(1,1,2) &
             + s(ij,1)*(THERMO_AI(1,1,1)-THERMO_AI(1,1,2))&
             + s(ij,2)*(THERMO_AI(1,1,3)-THERMO_AI(1,1,1))
        WMEAN_INV = WGHT_INV(2) + s(ij,1)*(WGHT_INV(1)-WGHT_INV(2)) - s(ij,2)*WGHT_INV(1)

        gama(ij) = HEAT_CAPACITY/(HEAT_CAPACITY-GRATIO*WMEAN_INV)

     ENDDO

! ###################################################################
! General mixture of NSP species. s contains NSP-1 mass fractions
! ###################################################################
  ELSE
     DO ij = 1, nx*ny*nz
! pass species to YMASS vector
        YMASS(NSP) = C_1_R
        DO is = 1, NSP-1
           YMASS(is) = s(ij,is)
           YMASS(NSP)= YMASS(NSP)-YMASS(is)
        ENDDO
! calculate cp and molecular weight of the mixture
        WMEAN_INV     = C_0_R
        HEAT_CAPACITY = C_0_R
        DO is = 1,NSP
           IF ( T(ij) .LT. THERMO_TLIM(3,is) ) THEN
              im = 2
           ELSE
              im = 1
           ENDIF
           HEAT_CAPACITY_I = C_0_R
           DO icp = NCP_CHEMKIN,1,-1
              HEAT_CAPACITY_I = HEAT_CAPACITY_I*T(ij) + THERMO_AI(icp,im,is)
           ENDDO
           HEAT_CAPACITY = HEAT_CAPACITY + YMASS(is)*HEAT_CAPACITY_I
           WMEAN_INV     = WMEAN_INV     + YMASS(is)*WGHT_INV(is)
        ENDDO

        gama(ij) = HEAT_CAPACITY/(HEAT_CAPACITY-GRATIO*WMEAN_INV)

     ENDDO

  ENDIF

  RETURN
END SUBROUTINE THERMO_GAMMA

