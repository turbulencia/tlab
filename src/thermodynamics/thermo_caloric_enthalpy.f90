#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/21 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Computing enthalpy from temperature and species mass fractions
!# according to the caloric equation of state.
!#
!# h = \sum Y_ih_i 
!#
!########################################################################
SUBROUTINE THERMO_CALORIC_ENTHALPY(nx,ny,nz, s, T, h)

  USE THERMO_VARS, ONLY : imixture
  USE THERMO_VARS, ONLY : NSP, NCP_CHEMKIN, THERMO_AI, THERMO_TLIM
  USE THERMO_VARS, ONLY : YMASS

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: T
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: H

! -------------------------------------------------------------------
  TINTEGER i, is, im, icp
  TREAL ENTHALPY_I, ENTHALPY_V, ENTHALPY_D, ENTHALPY_L

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO i = 1, nx*ny*nz
        h(i) = T(i)
     ENDDO

! ###################################################################
! Mixture defined by a conserved scalar Z, Y_i=f_i(Z)
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_BS          .OR.&
            imixture .EQ. MIXT_TYPE_BSZELDOVICH .OR.&
            imixture .EQ. MIXT_TYPE_QUASIBS    ) THEN
!     DO i = 1,nx*ny*nz
!
!#define MACRO_ZINPUT s(i,inb_scal)
!#include "dns_chem_mass.h"
!
!#define MACRO_TINPUT T(i)
!#include "dns_chem_enth.h"
!
!        h(i) = CH_H
!     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO i = 1, nx*ny*nz
        ENTHALPY_V = THERMO_AI(1,1,1)*T(i) + THERMO_AI(6,1,1)
        ENTHALPY_D = THERMO_AI(1,1,2)*T(i) + THERMO_AI(6,1,2)
        ENTHALPY_L = THERMO_AI(1,1,3)*T(i) + THERMO_AI(6,1,3)
        h(i) = (s(i,1)-s(i,2))*ENTHALPY_V + (C_1_R-s(i,1))*ENTHALPY_D &
             + s(i,2)*ENTHALPY_L
! only gas part
!        h(i) = ( (s(i,1)-s(i,2))*ENTHALPY_V + (C_1_R-s(i,1))*ENTHALPY_D )
! $           /(C_1_R-s(i,2))
     ENDDO

! ###################################################################
! General mixture of NSP species. s contains NSP-1 mass fractions
! ###################################################################
  ELSE
     h(1:nx*ny*nz) = C_0_R
     DO i = 1, nx*ny*nz
! pass species to YMASS vector
        YMASS(NSP) = C_1_R
        DO is = 1, NSP-1
           YMASS(is) = s(i,is)
           YMASS(NSP)= YMASS(NSP)-YMASS(is)
        ENDDO
! calculate enthalpy of the mixture
        DO is = 1,NSP
           IF ( T(i) .LT. THERMO_TLIM(3,is) ) THEN; im = 2
           ELSE;                                    im = 1; ENDIF
           ENTHALPY_I = C_0_R
           DO icp = NCP_CHEMKIN,1,-1
              ENTHALPY_I = ENTHALPY_I*T(i) + THERMO_AI(icp,im,is)/M_REAL(icp)
           ENDDO
           ENTHALPY_I = ENTHALPY_I*T(i) + THERMO_AI(6,im,is)
           h(i)       = h(i)      + YMASS(is)*ENTHALPY_I
        ENDDO
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE THERMO_CALORIC_ENTHALPY
