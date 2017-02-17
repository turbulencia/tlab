#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2007/11/08 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate qsat from rho and e
!#
!########################################################################
!# ARGUMENTS 
!#
!# qsat    In    contains first iteration on qsat
!#         Out   Value of qsat
!#
!########################################################################
SUBROUTINE THERMO_CALORIC_QSAT(nx,ny,nz, e,rho, qsat, T)

  USE THERMO_GLOBAL, ONLY : GRATIO, WGHT_INV, THERMO_AI, THERMO_PSAT, NPSAT

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL e(*), rho(*), qsat(*)
  TREAL T(*)

! -------------------------------------------------------------------
  TINTEGER i, is, inr, nrmax, ipsat
  TREAL LATENT_HEAT, HEAT_CAPACITY_VD, NEWTONRAPHSON_ERROR_LOC
  TREAL B_LOC(10), FUN, DER, B_LOC_CONST_2, B_LOC_CONST_3
  TREAL t_loc, HEAT_CAPACITY_DV, psat

! ###################################################################
  NEWTONRAPHSON_ERROR_LOC = 0
! maximum number of iterations in Newton-Raphson
  nrmax = 3

! reference case q_l = 0
  HEAT_CAPACITY_VD = THERMO_AI(1,1,1)-THERMO_AI(1,1,2)+GRATIO*WGHT_INV(2)-GRATIO*WGHT_INV(1)
  DO i = 1,nx*ny*nz
     T(i) = (e(i)-THERMO_AI(6,1,2)-qsat(i)*(THERMO_AI(6,1,1)-THERMO_AI(6,1,2)))/&
          (THERMO_AI(1,1,2)-GRATIO*WGHT_INV(2)+qsat(i)*HEAT_CAPACITY_VD)
  ENDDO

! initialize homogeneous data
  LATENT_HEAT      = THERMO_AI(6,1,1)-THERMO_AI(6,1,2)
  HEAT_CAPACITY_DV = THERMO_AI(1,1,2)-THERMO_AI(1,1,1)&
       + GRATIO*WGHT_INV(1) - GRATIO*WGHT_INV(2)
  DO i = 1,9
     B_LOC(i) =-THERMO_PSAT(i)*LATENT_HEAT
  ENDDO
  B_LOC(10) = C_0_R
  DO i = 2,10
     B_LOC(i) = B_LOC(i) + THERMO_PSAT(i-1)*HEAT_CAPACITY_DV
  ENDDO
  B_LOC_CONST_2 = B_LOC(2)
  B_LOC_CONST_3 = B_LOC(3)

! loop on all points
  DO i = 1, nx*ny*nz
     B_LOC(2) = B_LOC_CONST_2 + rho(i)*WGHT_INV(1)*&
          ( e(i)-THERMO_AI(6,1,2) )
     B_LOC(3) = B_LOC_CONST_3 - rho(i)*WGHT_INV(1)*&
          ( THERMO_AI(1,1,2)-GRATIO*WGHT_INV(2) )
! Newton-Raphson 
     t_loc = T(i)
     DO inr = 1,nrmax
        FUN = B_LOC(10)
        DER = C_0_R
        DO is = 9,1,-1
           FUN = FUN*t_loc + B_LOC(is)
           DER = DER*t_loc + B_LOC(is+1)*M_REAL(is)
        ENDDO
        t_loc = t_loc - FUN/DER
     ENDDO
     NEWTONRAPHSON_ERROR_LOC = MAX(NEWTONRAPHSON_ERROR_LOC,ABS(FUN/DER)/t_loc)
! calculate saturation specific humidity, in array z1(1,2).
! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid routine calls
     psat = C_0_R
     DO ipsat = NPSAT,1,-1
        psat = psat*t_loc + THERMO_PSAT(ipsat)
     ENDDO
     qsat(i) = psat/(rho(i)*t_loc*WGHT_INV(1))

  ENDDO

  RETURN
END SUBROUTINE THERMO_CALORIC_QSAT

