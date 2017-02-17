#include "types.h"
#include "dns_error.h"

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
!# Routint THERMO_POLYNOMIAL_PSAT is duplicated here to avoid array calls
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_PH(nx,ny,nz, s, p, h, T, dqldqt)

  USE DNS_CONSTANTS, ONLY : efile
  USE THERMO_GLOBAL, ONLY : MRATIO, WGHT_INV, THERMO_AI, THERMO_PSAT, NPSAT, dsmooth
  USE THERMO_GLOBAL, ONLY : NEWTONRAPHSON_ERROR
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                     INTENT(IN)   :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)   :: h,p
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(INOUT):: s        ! (*,1) is q_t, (*,2) is q_l
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT)  :: T
  TREAL  dqldqt(*)

! -------------------------------------------------------------------
  TINTEGER ij, is, inr, nrmax, ipsat
  TREAL psat, qsat, P_LOC, B_LOC(10), FUN, DER
  TREAL LATENT_HEAT, HEAT_CAPACITY_LD, HEAT_CAPACITY_LV, HEAT_CAPACITY_VD
  TREAL ALPHA_1, ALPHA_2, BETA_1, BETA_2, alpha, beta
  TREAL dummy

! ###################################################################
  IF ( dsmooth .GT. C_0_R ) THEN
     CALL IO_WRITE_ASCII(efile, 'THERMO_AIRWATER_PH. SmoothUndeveloped.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  NEWTONRAPHSON_ERROR = C_0_R
! maximum number of iterations in Newton-Raphson
  nrmax = 5

! initialize homogeneous data
  LATENT_HEAT      = THERMO_AI(6,1,1)-THERMO_AI(6,1,3)
  HEAT_CAPACITY_LV = THERMO_AI(1,1,3)-THERMO_AI(1,1,1)
  HEAT_CAPACITY_LD = THERMO_AI(1,1,3)-THERMO_AI(1,1,2)
  HEAT_CAPACITY_VD = HEAT_CAPACITY_LD - HEAT_CAPACITY_LV

  ALPHA_1 = WGHT_INV(2) /WGHT_INV(1) *LATENT_HEAT - THERMO_AI(6,1,2)
  ALPHA_2 = THERMO_AI(6,1,2) -THERMO_AI(6,1,1) &
          + LATENT_HEAT *( C_1_R -WGHT_INV(2) /WGHT_INV(1) )
  BETA_1  = WGHT_INV(2) /WGHT_INV(1) *HEAT_CAPACITY_LV + THERMO_AI(1,1,2)
  BETA_2  = HEAT_CAPACITY_LD - WGHT_INV(2) /WGHT_INV(1) *HEAT_CAPACITY_LV

! ###################################################################
  DO ij = 1, nx*ny*nz

     P_LOC = MRATIO*p(ij)
     
! -------------------------------------------------------------------
! reference case assuming q_l = 0
! -------------------------------------------------------------------
     T(ij) = (h(ij)-THERMO_AI(6,1,2)-s(ij,1)*(THERMO_AI(6,1,1)-THERMO_AI(6,1,2)))/&
          (THERMO_AI(1,1,2)+s(ij,1)*(THERMO_AI(1,1,1)-THERMO_AI(1,1,2)))

! calculate saturation specific humidity q_sat(T,p) in array s(1,2)
     psat = C_0_R
     DO ipsat = NPSAT,1,-1
        psat = psat*T(ij) + THERMO_PSAT(ipsat)
     ENDDO
     dummy = WGHT_INV(2) /WGHT_INV(1) /( P_LOC/psat -C_1_R )
     qsat = dummy /( C_1_R +dummy )

     IF ( qsat .GE. s(ij,1) ) THEN
        s(ij,2) = C_0_R

! -------------------------------------------------------------------
! if q_s < q_t, then we have to recalculate T
! -------------------------------------------------------------------
     ELSE
! preparing Newton-Raphson 
        alpha = (ALPHA_1 + s(ij,1)*ALPHA_2 + h(ij)) /P_LOC
        beta  = (BETA_1  + s(ij,1)*BETA_2         ) /P_LOC
        B_LOC(1) = h(ij)-THERMO_AI(6,1,2) - s(ij,1)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,2)) -&
             THERMO_PSAT(1)*alpha
        DO is = 2,9
           B_LOC(is) = THERMO_PSAT(is-1)*beta - THERMO_PSAT(is)*alpha
        ENDDO
        B_LOC(2)  = B_LOC(2) - THERMO_AI(1,1,2) - s(ij,1)*HEAT_CAPACITY_LD
        B_LOC(10) = THERMO_PSAT(9)*beta

! executing Newton-Raphson 
        DO inr = 1,nrmax
           FUN = B_LOC(10)
           DER = C_0_R
           DO is = 9,1,-1
              FUN = FUN*T(ij) + B_LOC(is)
              DER = DER*T(ij) + B_LOC(is+1)*M_REAL(is)
           ENDDO
           T(ij) = T(ij) - FUN/DER
        ENDDO
        NEWTONRAPHSON_ERROR = MAX(NEWTONRAPHSON_ERROR,ABS(FUN/DER)/T(ij))

! recalculate saturation specific humidity q_vs(T,p,qt)
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*T(ij) + THERMO_PSAT(ipsat)
        ENDDO
        dummy = WGHT_INV(2) /WGHT_INV(1) /( P_LOC/psat -C_1_R )
        s(ij,2) = dummy *( C_1_R -s(ij,1) )

! liquid content
        s(ij,2) = s(ij,1)-s(ij,2)

     ENDIF

  ENDDO


#ifdef USE_MPI
  CALL MPI_ALLREDUCE&
       (NEWTONRAPHSON_ERROR, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  NEWTONRAPHSON_ERROR = dummy
#endif

  RETURN
END SUBROUTINE THERMO_AIRWATER_PH

!########################################################################
!########################################################################
! iterative method based on (rho,e,q_i)
SUBROUTINE THERMO_AIRWATER_PH2(nx, ny, nz, z1, p, h, T)

  USE THERMO_GLOBAL, ONLY : GRATIO, MRATIO, THERMO_AI

  IMPLICIT NONE

#include "integers.h"
  
  TINTEGER nx, ny, nz
  TREAL z1(nx*ny*nz,*), T(*), h(*), p(*)

! -------------------------------------------------------------------
  TINTEGER ij, iter, niter
  TREAL r_loc, e_loc, t_loc, z1_loc(2), dummy, prefactor

! ###################################################################
  niter = 5
  prefactor = GRATIO*MRATIO ! = (gama0-C_1_R)*mach*mach

  DO ij = 1,nx*ny*nz
! -------------------------------------------------------------------
! initialize, q_l=0
! -------------------------------------------------------------------
     z1_loc(1) = z1(ij,1)
     z1_loc(2) = C_0_R
     t_loc = (h(ij)-THERMO_AI(6,1,2)-z1(ij,1)*(THERMO_AI(6,1,1)-THERMO_AI(6,1,2)))/&
          (THERMO_AI(1,1,2)+z1(ij,1)*(THERMO_AI(1,1,1)-THERMO_AI(1,1,2)))

! -------------------------------------------------------------------
! iteration
! -------------------------------------------------------------------
     DO iter = 1,niter
! calculate density from temperature/composition
        CALL THERMO_THERMAL_DENSITY(i1, i1, i1, z1_loc, p(ij), t_loc, r_loc)

! calculate energy
        e_loc = h(ij) - prefactor*p(ij)/r_loc

! solve equilibrium (rho,e,q_i)
        CALL THERMO_AIRWATER_RE(i1, i1, i1, z1_loc, e_loc, r_loc, t_loc, dummy)

     ENDDO
     z1(ij,2) = z1_loc(2)
     T(ij)    = t_loc

  ENDDO

  RETURN
END SUBROUTINE THERMO_AIRWATER_PH2
