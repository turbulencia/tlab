!########################################################################
!# Tool/Library DNS
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
!# Computing energy from temperature and species mass fractions
!# according to the caloric equation of state.
!#
!# e=\sum Y_ih_i - TR_0\sum Y_i/W_i
!#
!# which in nondimensional form is
!#
!# e=\sum Y_ih_i - \sum Y_i T*(R_0/C_{p,0]W_0)*W_0/W_i
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

SUBROUTINE THERMO_CALORIC_ENERGY(nx,ny,nz, s,T, e)

  USE THERMO_GLOBAL, ONLY : imixture, gama0, GRATIO
  USE THERMO_GLOBAL, ONLY : NSP, NCP_CHEMKIN, WGHT_INV, THERMO_AI, THERMO_TLIM
  USE THERMO_GLOBAL, ONLY : YMASS

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: T
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: e

! -------------------------------------------------------------------
  TINTEGER i, is, im, icp
  TREAL ENTHALPY_I, WMEAN_INV
  TREAL ENERGY_V, ENERGY_D, ENERGY_L

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO i = 1, nx*ny*nz
        e(i) = T(i)/gama0
     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO i = 1, nx*ny*nz
        ENERGY_V = (THERMO_AI(1,1,1)- GRATIO*WGHT_INV(1))*T(i) + THERMO_AI(6,1,1)
        ENERGY_D = (THERMO_AI(1,1,2)- GRATIO*WGHT_INV(2))*T(i) + THERMO_AI(6,1,2)
        ENERGY_L = THERMO_AI(1,1,3)*T(i) + THERMO_AI(6,1,3)
        e(i) = (s(i,1)-s(i,2))*ENERGY_V + (C_1_R-s(i,1))*ENERGY_D + s(i,2)*ENERGY_L
     ENDDO

! ###################################################################
! General mixture of NSP species. s contains NSP-1 mass fractions
! ###################################################################
  ELSE
     e(1:nx*ny*nz) = C_0_R
     DO i = 1, nx*ny*nz
! pass species to YMASS vector
        YMASS(NSP) = C_1_R
        DO is = 1, NSP-1
           YMASS(is) = s(i,is)
           YMASS(NSP)= YMASS(NSP)-YMASS(is)
        ENDDO
! calculate enthalpy of the mixture
        WMEAN_INV = C_0_R
        DO is = 1,NSP
           IF ( T(i) .LT. THERMO_TLIM(3,is) ) THEN
              im = 2
           ELSE
              im = 1
           ENDIF
           ENTHALPY_I = C_0_R
           DO icp = NCP_CHEMKIN,1,-1
              ENTHALPY_I = ENTHALPY_I*T(i) + THERMO_AI(icp,im,is)/M_REAL(icp)
           ENDDO
           ENTHALPY_I = ENTHALPY_I*T(i) + THERMO_AI(6,im,is)
           e(i)       = e(i)      + YMASS(is)*ENTHALPY_I
           WMEAN_INV  = WMEAN_INV + YMASS(is)*WGHT_INV(is)
        ENDDO
! go from enthalpy to energy
        e(i) = e(i) - GRATIO*WMEAN_INV*T(i)
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE THERMO_CALORIC_ENERGY
