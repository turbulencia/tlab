#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/11/08 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Extracted from THERMO_AIRWATER_RE. Case dsmooth > 0 not yet validated !
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_RH(nx, ny, nz, z1, h, rho, T, dqldqt)

  USE THERMO_GLOBAL, ONLY : GRATIO, WGHT_INV, THERMO_AI, THERMO_PSAT, NPSAT, dsmooth
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx, ny, nz
  TREAL h(*), z1(nx*ny*nz,*), rho(*), dqldqt(*)
  TREAL T(*)

! -------------------------------------------------------------------
  TINTEGER i, is, inr, nrmax, ipsat
  TREAL LATENT_HEAT, HEAT_CAPACITY_LD, HEAT_CAPACITY_LV, HEAT_CAPACITY_VD
  TREAL B_LOC(10), FUN, DER, B_LOC_CONST_2, B_LOC_CONST_3
  TREAL dsmooth_loc
  TREAL t_loc, qsat, HEAT_CAPACITY_DV, psat, ERROR_LOC
  TREAL alpha, dummy1, dummy2
#ifdef USE_MPI
  TREAL dummy
#endif

! ###################################################################
  ERROR_LOC = C_0_R
! maximum number of iterations in Newton-Raphson
  nrmax = 3

! reference case q_l = 0
  DO i = 1,nx*ny*nz
     T(i) = (h(i)-THERMO_AI(6,1,2)-z1(i,1)*(THERMO_AI(6,1,1)-THERMO_AI(6,1,2)))/&
          (THERMO_AI(1,1,2)+z1(i,1)*(THERMO_AI(1,1,1)-THERMO_AI(1,1,2)))
  ENDDO

! -------------------------------------------------------------------
! calculate saturation specific humidity q_s(\rho,e)
! -------------------------------------------------------------------
  IF ( dsmooth .LE. C_0_R ) THEN
! calculate saturation specific humidity, in array z1(1,2).
! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid array calls
     DO i = 1,nx*ny*nz
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*T(i) + THERMO_PSAT(ipsat)
        ENDDO
        z1(i,2) = psat/(rho(i)*T(i)*WGHT_INV(1))
     ENDDO

  ELSE
! initialize homogeneous data
     LATENT_HEAT      = THERMO_AI(6,1,1)-THERMO_AI(6,1,2)
     HEAT_CAPACITY_DV = THERMO_AI(1,1,2)-THERMO_AI(1,1,1)
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
             ( h(i)-THERMO_AI(6,1,2) )
        B_LOC(3) = B_LOC_CONST_3 - rho(i)*WGHT_INV(1)*&
             ( THERMO_AI(1,1,2) )
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
        ERROR_LOC = MAX(ERROR_LOC,ABS(FUN/DER)/t_loc)
! calculate saturation specific humidity, in array z1(1,2).
! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid routine calls
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*t_loc + THERMO_PSAT(ipsat)
        ENDDO
        z1(i,2) = psat/(rho(i)*t_loc*WGHT_INV(1))

! calculate dqldqt
        qsat = z1(i,2)
        alpha = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) -&
             (THERMO_AI(1,1,3)-THERMO_AI(1,1,1))*t_loc
        alpha = alpha/(t_loc*GRATIO*WGHT_INV(1)) - C_1_R

        dummy1 = THERMO_AI(6,1,2) - THERMO_AI(6,1,3) -&
             (THERMO_AI(1,1,3)-THERMO_AI(1,1,2))*t_loc
        dummy1 = dummy1*qsat*alpha

        dummy2 = THERMO_AI(6,1,2) - THERMO_AI(6,1,1) +&
             t_loc*GRATIO*WGHT_INV(1)*alpha*(alpha+C_1_R)
        dummy2 = dummy2*qsat+h(i)-THERMO_AI(6,1,2)

        dqldqt(i) = C_1_R-dummy1/dummy2
     ENDDO

  ENDIF

! -------------------------------------------------------------------
! calculate final T and q_l
! -------------------------------------------------------------------
! initialize homogeneous data
  LATENT_HEAT      = THERMO_AI(6,1,1)-THERMO_AI(6,1,3)
  HEAT_CAPACITY_LV = THERMO_AI(1,1,3)-THERMO_AI(1,1,1)
  HEAT_CAPACITY_LD = THERMO_AI(1,1,3)-THERMO_AI(1,1,2)
  HEAT_CAPACITY_VD = HEAT_CAPACITY_LD - HEAT_CAPACITY_LV
  DO i = 1,9
     B_LOC(i) =-THERMO_PSAT(i)*LATENT_HEAT
  ENDDO
  B_LOC(10) = C_0_R
  DO i = 2,10
     B_LOC(i) = B_LOC(i) + THERMO_PSAT(i-1)*HEAT_CAPACITY_LV
  ENDDO
  B_LOC_CONST_2 = B_LOC(2)
  B_LOC_CONST_3 = B_LOC(3)

! loop on all points
  DO i = 1, nx*ny*nz
     qsat = z1(i,2)

     IF ( qsat .GE. z1(i,1) ) THEN
        z1(i,2) = C_0_R
        IF ( dsmooth .GT. C_0_R ) THEN
           dsmooth_loc = dsmooth*qsat
           z1(i,2) = dsmooth_loc*dqldqt(i)&
                *LOG(EXP((z1(i,1)-qsat)/dsmooth_loc)+C_1_R)
! change T consistently
           T(i) = ( h(i) - &
                z1(i,1)*(THERMO_AI(6,1,1)-THERMO_AI(6,1,2)) - THERMO_AI(6,1,2) +&
                z1(i,2)*LATENT_HEAT )/&
                ( z1(i,1)*HEAT_CAPACITY_VD + THERMO_AI(1,1,2)&
                + z1(i,2)*HEAT_CAPACITY_LV )
        ENDIF

! if q_s < q_t, then we have to repeat calculation of T
     ELSE
        B_LOC(2) = B_LOC_CONST_2 + rho(i)*WGHT_INV(1)*&
             ( h(i)-z1(i,1)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,2))-THERMO_AI(6,1,2) )
        B_LOC(3) = B_LOC_CONST_3 - rho(i)*WGHT_INV(1)*&
             ( z1(i,1)*HEAT_CAPACITY_LD + THERMO_AI(1,1,2) )
! IF ( dsmooth .GT. C_0_R ) THEN
! dsmooth_loc = dsmooth*qsat
! alpha =-( dsmooth_loc*LOG(EXP((z1(i,1)-qsat)/dsmooth_loc)+C_1_R) 
! $                 -(z1(i,1)-qsat) )*dqldqt(i) 
! B_LOC(2) = B_LOC(2) - rho(i)*WGHT_INV(1)*alpha*LATENT_HEAT
! B_LOC(3) = B_LOC(3) + rho(i)*WGHT_INV(1)*alpha*HEAT_CAPACITY_LV
! ENDIF

! Newton-Raphson 
        DO inr = 1,nrmax
           FUN = B_LOC(10)
           DER = C_0_R
           DO is = 9,1,-1
              FUN = FUN*T(i) + B_LOC(is)
              DER = DER*T(i) + B_LOC(is+1)*M_REAL(is)
           ENDDO
           T(i) = T(i) - FUN/DER
        ENDDO
        ERROR_LOC = MAX(ERROR_LOC,ABS(FUN/DER)/T(i))

! calculate saturation specific humidity, in array z1(1,2).
! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid routine calls
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*T(i) + THERMO_PSAT(ipsat)
        ENDDO
        z1(i,2) = psat/(rho(i)*T(i)*WGHT_INV(1))

! liquid content
        z1(i,2) = z1(i,1)-z1(i,2)
        IF ( dsmooth .GT. C_0_R ) THEN
           dsmooth_loc = dsmooth*qsat
           z1(i,2) = dsmooth_loc*dqldqt(i)&
                *LOG(EXP((z1(i,1)-qsat)/dsmooth_loc)+C_1_R) &
                + z1(i,2) - dqldqt(i)*(z1(i,1)-qsat)
! change T consistently
           T(i) = ( h(i) - &
                z1(i,1)*(THERMO_AI(6,1,1)-THERMO_AI(6,1,2)) - THERMO_AI(6,1,2) +&
                z1(i,2)*LATENT_HEAT )/&
                ( z1(i,1)*HEAT_CAPACITY_VD + THERMO_AI(1,1,2)&
                + z1(i,2)*HEAT_CAPACITY_LV )
        ENDIF
     ENDIF

  ENDDO

#ifdef USE_MPI
  CALL MPI_ALLREDUCE(ERROR_LOC, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  ERROR_LOC = dummy
#endif

  RETURN
END SUBROUTINE THERMO_AIRWATER_RH
