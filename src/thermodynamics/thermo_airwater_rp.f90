#include "types.h"

!########################################################################
!# Tool DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/08/02 - J.P. Mellado
!#              Created
!# 2007/09/24 - J.P. Mellado
!#              Adding smoothing factor
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate T and liquid content from rho, p and water content 
!# using thermal equation of state.
!# The difference with THERMO_THERMAL_TEMPERATURE is that the equilibrium
!# partition of q_t between vapor and liquid here is not known and q_l
!# has to be calculated.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_RP(nx, ny, nz, z1, p, rho, T, dqldqt)

  USE THERMO_VARS, ONLY : GRATIO, MRATIO, WGHT_INV, THERMO_AI, THERMO_PSAT, NPSAT, dsmooth

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL z1(nx*ny*nz,*), T(*), rho(*), p(*), dqldqt(*)

! -------------------------------------------------------------------
  TINTEGER ij, is, inr, nrmax, ipsat
  TREAL psat, qsat, B_LOC(9), FUN, DER, ERROR_LOC, dsmooth_loc
  TREAL t_loc, B_LOC_CONST_1, B_LOC_CONST_2, dummy, alpha

! ###################################################################
  ERROR_LOC = C_0_R
! maximum number of iterations in Newton-Raphson
  nrmax = 3

! reference case q_l = 0
  DO ij = 1,nx*ny*nz
     T(ij) = MRATIO*p(ij)/&
          ( (WGHT_INV(2)+z1(ij,1)*(WGHT_INV(1)-WGHT_INV(2))) * rho(ij) )
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
        z1(ij,2) = psat/(rho(ij)*T(ij)*WGHT_INV(1))
     ENDDO

  ELSE
! initialize homogeneous data
     DO ij = 1,9
        B_LOC(ij) = THERMO_PSAT(ij)*(WGHT_INV(2)/WGHT_INV(1)-C_1_R)
     ENDDO
     B_LOC_CONST_1 = B_LOC(1)
     B_LOC_CONST_2 = B_LOC(2)

! loop on all points
     DO ij = 1, nx*ny*nz
        B_LOC(1) = B_LOC_CONST_1 + MRATIO*p(ij)
        B_LOC(2) = B_LOC_CONST_2 - rho(ij)*WGHT_INV(2)

! Newton-Raphson 
        t_loc = T(ij)
        DO inr = 1,nrmax
           FUN = B_LOC(9)
           DER = C_0_R
           DO is = 8,1,-1
              FUN = FUN*t_loc + B_LOC(is)
              DER = DER*t_loc + B_LOC(is+1)*M_REAL(is)
           ENDDO
           t_loc = t_loc - FUN/DER
        ENDDO
        ERROR_LOC = MAX(ERROR_LOC,ABS(FUN/DER)/t_loc)

! calculate saturation specific humidity, in array z1(1,2).
        z1(ij,2) = (MRATIO*p(ij)-rho(ij)*WGHT_INV(2)*t_loc)/&
             (rho(ij)*(WGHT_INV(1)-WGHT_INV(2))*t_loc)

! calculate dqldqt
        qsat = z1(ij,2)
        alpha = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) -&
             (THERMO_AI(1,1,3)-THERMO_AI(1,1,1)+GRATIO*WGHT_INV(1))*t_loc
        alpha = alpha/(t_loc*GRATIO*WGHT_INV(1)) - C_1_R

        dummy = MRATIO*p(ij)/(qsat*rho(ij)*WGHT_INV(1)*t_loc)

        dqldqt(ij) = C_1_R - alpha*WGHT_INV(2)/WGHT_INV(1)/(C_1_R+dummy)
     ENDDO

  ENDIF

! -------------------------------------------------------------------
! check if condesation occurs and, if so, solve nonlinear equation
! -------------------------------------------------------------------
! initialize homogeneous data
  DO ij = 1,9
     B_LOC(ij) = THERMO_PSAT(ij)
  ENDDO

! loop on all points
  DO ij = 1, nx*ny*nz
     qsat = z1(ij,2)

     IF ( qsat .GT. z1(ij,1) ) THEN
        z1(ij,2) = C_0_R
        IF ( dsmooth .GT. C_0_R ) THEN
           dsmooth_loc = dsmooth*qsat
           z1(ij,2) = dsmooth_loc*dqldqt(ij)&
                *LOG(EXP((z1(ij,1)-qsat)/dsmooth_loc)+C_1_R)
! change T consistently
           T(ij) = MRATIO*p(ij)/&
                ( ((z1(ij,1)-z1(ij,2))*WGHT_INV(1)+(C_1_R-z1(ij,1))*WGHT_INV(2)) * rho(ij) )
        ENDIF

! if q_s < q_t, then we have to repeat calculation of T
     ELSE
        B_LOC(1) = THERMO_PSAT(1) - MRATIO*p(ij)
        B_LOC(2) = THERMO_PSAT(2) + (C_1_R-z1(ij,1))*rho(ij)*WGHT_INV(2)

! Newton-Raphson 
        DO inr = 1,nrmax
           FUN = B_LOC(9)
           DER = C_0_R
           DO is = 8,1,-1
              FUN = FUN*T(ij) + B_LOC(is)
              DER = DER*T(ij) + B_LOC(is+1)*M_REAL(is)
           ENDDO
           T(ij) = T(ij) - FUN/DER
        ENDDO
        ERROR_LOC = MAX(ERROR_LOC,ABS(FUN/DER)/T(ij))

! calculate saturation specific humidity, in array z1(1,2).
! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid array calls
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*T(ij) + THERMO_PSAT(ipsat)
        ENDDO
        z1(ij,2) = psat/(rho(ij)*T(ij)*WGHT_INV(1))

! liquid content
        z1(ij,2) = z1(ij,1) - z1(ij,2)
        IF ( dsmooth .GT. C_0_R ) THEN
           dsmooth_loc = dsmooth*qsat
           z1(ij,2) = dsmooth_loc*dqldqt(ij)&
                *LOG(EXP((z1(ij,1)-qsat)/dsmooth_loc)+C_1_R) &
                + z1(ij,2) - dqldqt(ij)*(z1(ij,1)-qsat)
! change T consistently
           T(ij) = MRATIO*p(ij)/&
                ( ((z1(ij,1)-z1(ij,2))*WGHT_INV(1)+(C_1_R-z1(ij,1))*WGHT_INV(2)) * rho(ij) )
        ENDIF
     ENDIF

  ENDDO

  RETURN
END SUBROUTINE THERMO_AIRWATER_RP
