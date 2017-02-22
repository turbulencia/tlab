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
!# Iterative method based on (rho,e,q_i)
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_PH_RE(nx, ny, nz, z1, p, h, T)

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
END SUBROUTINE THERMO_AIRWATER_PH_RE
