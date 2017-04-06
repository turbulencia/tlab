#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2007/04/07 - Alberto de Lozar
!#              Created
!# 2016/10/15 - J.P. Mellado
!#              Cleaning
!# 2017/02/17 - J.P. Mellado
!#              Adapting for background p profile
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculating the equilibrium T and q_l for given enthalpy and pressure.
!# Assumes often that THERMO_AI(6,1,1) = THERMO_AI(6,1,2) = 0
!#
!# Routine THERMO_POLYNOMIAL_PSAT is duplicated here to avoid array calls
!#
!# Smoothing according to Eq. 25 in Mellado et al., TCFD, 2010
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_PH(nx,ny,nz, s,h, e,p)

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
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(INOUT):: s        ! (*,1) is q_t, (*,2) is q_l
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)   :: h
  TREAL, DIMENSION(*),          INTENT(IN)   :: e,p

! -------------------------------------------------------------------
  TINTEGER ij, is, inr, nrmax, ipsat, i, jk
  TREAL psat, qsat, T_LOC, E_LOC, P_LOC, B_LOC(10), FUN, DER
  TREAL LATENT_HEAT, HEAT_CAPACITY_LD, HEAT_CAPACITY_LV, HEAT_CAPACITY_VD
  TREAL ALPHA_1, ALPHA_2, BETA_1, BETA_2, alpha, beta, rd_ov_rv
  TREAL dummy
  TREAL qvequ, dsmooth_loc, dqldqt, dqsdt
  
  TREAL Cd, Cdv, Lv0, Cvl

! ###################################################################
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

  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  rd_ov_rv = WGHT_INV(2) /WGHT_INV(1)

! ###################################################################
  ij = 0
  DO jk = 0,ny*nz-1
     is = MOD(jk,ny) +1
     P_LOC = MRATIO*p(is)
     E_LOC = e(is)
     
     DO i = 1,nx
        ij = ij +1
        
! -------------------------------------------------------------------
! reference case assuming ql = 0
! -------------------------------------------------------------------
        s(ij,2) = C_0_R
        T_LOC   = ( h(ij) -E_LOC ) /( Cd +s(ij,1) *Cdv )
        
! calculate saturation specific humidity q_sat(T,p)
        psat = C_0_R
        DO ipsat = NPSAT,1,-1
           psat = psat*T_LOC + THERMO_PSAT(ipsat)
        ENDDO
        dummy = rd_ov_rv /( P_LOC/psat -C_1_R )
        qsat = dummy /( C_1_R +dummy )
     
! -------------------------------------------------------------------
! calculate smoothed piecewise linear contribution, if needed
! -------------------------------------------------------------------
        IF ( dsmooth .GT. C_0_R ) THEN
! calculate dqsdt from dpsatdt (qs here mean qvequ)
           dqsdt = C_0_R
           DO ipsat = NPSAT-1,1,-1
              dqsdt = dqsdt *T_LOC + THERMO_PSAT(ipsat+1) *M_REAL(ipsat)
           ENDDO
           dqsdt = qsat / psat /( C_1_R - psat/P_LOC ) *dqsdt
! calculate dqldqt
           dqsdt = dqsdt / ( Cd + qsat *Cdv )
           dqldqt = ( C_1_R/( C_1_R - qsat ) + Cdv *T_LOC *dqsdt )/ &
                    ( C_1_R + ( Lv0 - Cvl *T_LOC ) *dqsdt )
           
           dsmooth_loc  = dsmooth *qsat
           s(ij,2) = dqldqt *dsmooth_loc *LOG( EXP( (s(ij,1)-qsat) /dsmooth_loc ) +C_1_R )

        ENDIF
        
! -------------------------------------------------------------------
! if q_s < q_t, then we have to recalculate T
! -------------------------------------------------------------------
        IF ( qsat .LT. s(ij,1) ) THEN
! preparing Newton-Raphson 
           alpha = (ALPHA_1 + s(ij,1)*ALPHA_2 + h(ij) - E_LOC) /P_LOC
           beta  = (BETA_1  + s(ij,1)*BETA_2                 ) /P_LOC
           B_LOC(1) = h(ij) -E_LOC + s(ij,1) *Lv0 - THERMO_PSAT(1)*alpha
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
                 FUN = FUN*T_LOC + B_LOC(is)
                 DER = DER*T_LOC + B_LOC(is+1)*M_REAL(is)
              ENDDO
              T_LOC = T_LOC - FUN/DER
           ENDDO
           NEWTONRAPHSON_ERROR = MAX(NEWTONRAPHSON_ERROR,ABS(FUN/DER)/T_LOC)
           
! calculate equilibrium vapor specific humidity
           psat = C_0_R
           DO ipsat = NPSAT,1,-1
              psat = psat*T_LOC + THERMO_PSAT(ipsat)
           ENDDO
           dummy = rd_ov_rv /( P_LOC/psat -C_1_R )
           qvequ = dummy *( C_1_R -s(ij,1) )

           IF ( dsmooth .GT. C_0_R ) THEN ! add correction
              s(ij,2) = s(ij,2) +s(ij,1) -qvequ - (s(ij,1) -qsat) *dqldqt
           ELSE                           ! or calculate new
              s(ij,2) =          s(ij,1) -qvequ
           ENDIF
           
        ENDIF

     ENDDO
  ENDDO

#ifdef USE_MPI
  CALL MPI_ALLREDUCE&
       (NEWTONRAPHSON_ERROR, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  NEWTONRAPHSON_ERROR = dummy
#endif

  RETURN
END SUBROUTINE THERMO_AIRWATER_PH
