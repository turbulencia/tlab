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

! ###################################################################
! ###################################################################
SUBROUTINE THERMO_ANELASTIC_AIRWATER_PH_RE(nx,ny,nz, s, e,p, wrk3d)

  USE THERMO_GLOBAL, ONLY : GRATIO, MRATIO

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(INOUT)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: wrk3d

! -------------------------------------------------------------------
  TINTEGER ij, jk, is, i, iter, niter
  TREAL r_loc, en_loc, t_loc, z1_loc(2), dummy, prefactor
  TREAL p_loc, e_loc

! ###################################################################
  niter = 5
  prefactor = GRATIO *MRATIO ! = (gama0-C_1_R)*mach*mach

  print*,'hello'
  
  s(:,3) = C_0_R ! initialize, q_l=0

  DO iter = 1,niter ! iteration
! calculate density in wrk3d
     CALL THERMO_ANELASTIC_DENSITY(nx,ny,nz, s, e,p, wrk3d)

! calculate energy
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)

         DO i = 1,nx
           ij = ij +1

           r_loc = wrk3d(ij)
           z1_loc(1) = s(ij,2)
           z1_loc(2) = s(ij,3)
           en_loc = s(ij,1) - E_LOC - prefactor *P_LOC /r_loc
           
! solve equilibrium (rho,e,q_i)
           CALL THERMO_AIRWATER_RE(i1,i1,i1, z1_loc, en_loc, r_loc, t_loc, dummy)
           
           s(ij,3) = z1_loc(2)
           
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_AIRWATER_PH_RE
