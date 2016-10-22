#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/10/05 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculating the equilibrium T and q_l for given enthalpy and pressure.
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_PH(nx, ny, nz, z1, p, h, T, dqldqt)

  USE DNS_CONSTANTS, ONLY : efile
  USE THERMO_GLOBAL, ONLY : MRATIO, WGHT_INV, THERMO_AI, THERMO_PSAT, NPSAT, dsmooth

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL z1(nx*ny*nz,*), T(*), h(*), p(*), dqldqt(*)

! -------------------------------------------------------------------
  TINTEGER ij, is, inr, nrmax, ipsat
  TREAL psat, qsat, B_LOC(10), FUN, DER, ERROR_LOC
  TREAL LATENT_HEAT, HEAT_CAPACITY_LD, HEAT_CAPACITY_LV, HEAT_CAPACITY_VD
  TREAL ALPHA_1, ALPHA_2, BETA_1, BETA_2, alpha, beta
  TREAL dummy

! ###################################################################
  ERROR_LOC = C_0_R
! maximum number of iterations in Newton-Raphson
  nrmax = 5

! reference case q_l = 0
  DO ij = 1,nx*ny*nz
     T(ij) = (h(ij)-THERMO_AI(6,1,2)-z1(ij,1)*(THERMO_AI(6,1,1)-THERMO_AI(6,1,2)))/&
          (THERMO_AI(1,1,2)+z1(ij,1)*(THERMO_AI(1,1,1)-THERMO_AI(1,1,2)))
  ENDDO

! -------------------------------------------------------------------
! calculate saturation specific humidity q_s(\rho,p)
! -------------------------------------------------------------------
  IF ( dsmooth .LE. C_0_R ) THEN
! calculate saturation specific humidity, in array z1(1,2).
! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid array calls
     DO ij = 1,nx*ny*nz
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*T(ij) + THERMO_PSAT(ipsat)
        ENDDO
        dummy = C_1_R/(MRATIO*p(ij)/psat-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        z1(ij,2) = dummy/(C_1_R+dummy)
     ENDDO
  ELSE
     CALL IO_WRITE_ASCII(efile, 'THERMO_AIRWATER_PH. SmoothUndeveloped.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! -------------------------------------------------------------------
! calculate final T and q_l
! -------------------------------------------------------------------
! initialize homogeneous data
  LATENT_HEAT      = THERMO_AI(6,1,1)-THERMO_AI(6,1,3)
  HEAT_CAPACITY_LV = THERMO_AI(1,1,3)-THERMO_AI(1,1,1)
  HEAT_CAPACITY_LD = THERMO_AI(1,1,3)-THERMO_AI(1,1,2)
  HEAT_CAPACITY_VD = HEAT_CAPACITY_LD - HEAT_CAPACITY_LV

  ALPHA_1 = WGHT_INV(2)/WGHT_INV(1)*LATENT_HEAT - THERMO_AI(6,1,2)
  ALPHA_2 = THERMO_AI(6,1,2)-THERMO_AI(6,1,1) +&
       LATENT_HEAT*(C_1_R-WGHT_INV(2)/WGHT_INV(1))
  BETA_1  = WGHT_INV(2)/WGHT_INV(1)*HEAT_CAPACITY_LV + THERMO_AI(1,1,2)
  BETA_2  = HEAT_CAPACITY_LD - WGHT_INV(2)/WGHT_INV(1)*HEAT_CAPACITY_LV

! loop on all points
  DO ij = 1, nx*ny*nz
     qsat = z1(ij,2)

     IF ( qsat .GE. z1(ij,1) ) THEN
        z1(ij,2) = C_0_R

! if q_s < q_t, then we have to repeat calculation of T
     ELSE
        alpha = (ALPHA_1 + z1(ij,1)*ALPHA_2 + h(ij))/(MRATIO*p(ij))
        beta  = (BETA_1  + z1(ij,1)*BETA_2         )/(MRATIO*p(ij))
        B_LOC(1) = h(ij)-THERMO_AI(6,1,2) - z1(ij,1)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,2)) -&
             THERMO_PSAT(1)*alpha
        DO is = 2,9
           B_LOC(is) = THERMO_PSAT(is-1)*beta - THERMO_PSAT(is)*alpha
        ENDDO
        B_LOC(2)  = B_LOC(2) - THERMO_AI(1,1,2) - z1(ij,1)*HEAT_CAPACITY_LD
        B_LOC(10) = THERMO_PSAT(9)*beta

! Newton-Raphson 
        DO inr = 1,nrmax
           FUN = B_LOC(10)
           DER = C_0_R
           DO is = 9,1,-1
              FUN = FUN*T(ij) + B_LOC(is)
              DER = DER*T(ij) + B_LOC(is+1)*M_REAL(is)
           ENDDO
           T(ij) = T(ij) - FUN/DER
        ENDDO
        ERROR_LOC = MAX(ERROR_LOC,ABS(FUN/DER)/T(ij))

! calculate saturation specific humidity, in array z1(1,2).
! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid routine calls
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*T(ij) + THERMO_PSAT(ipsat)
        ENDDO
        dummy = C_1_R/(MRATIO*p(ij)/psat-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        z1(ij,2) = dummy*(C_1_R-z1(ij,1))

! liquid content
        z1(ij,2) = z1(ij,1)-z1(ij,2)

     ENDIF

  ENDDO

  RETURN
END SUBROUTINE THERMO_AIRWATER_PH
