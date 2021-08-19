#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/09 - J.P. Mellado
!#              Created
!# 2007/06/21 - J.P. Mellado
!#              Introducing case AIRWATER
!#
!########################################################################
!# DESCRIPTION
!#
!# Computing temperature from density, energy and species mass fractions:
!#
!# e=\sum Y_ih_i - TR_0\sum Y_i/W_i
!#
!# I need rho for the case in which the equilibrium composition is
!# to be calculated.
!#
!########################################################################
SUBROUTINE THERMO_CALORIC_TEMPERATURE(nx,ny,nz, s,e,rho, T, wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_PROCS
  USE THERMO_VARS, ONLY : imixture, gama0, GRATIO
  USE THERMO_VARS, ONLY : NSP, NCP_CHEMKIN, WGHT_INV, THERMO_AI

  IMPlICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)    :: e, rho
  TREAL, DIMENSION(nx*ny*nz),   INTENT(INOUT) :: wrk3d
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT)   :: T

! -------------------------------------------------------------------
  TINTEGER i, is
  TREAL WMEAN_INV, HEAT_CAPACITY, FORMATION_ENERGY

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO i = 1, nx*ny*nz
        T(i) = gama0*e(i)
     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     CALL THERMO_AIRWATER_RE(nx, ny, nz, s, e, rho, T, wrk3d)

! ###################################################################
! General mixture of NSP species. s contains NSP-1 mass fractions
! ###################################################################
  ELSE
! -------------------------------------------------------------------
! Cp linear with T. It assumes also only one temperature range
! -------------------------------------------------------------------
     IF ( NCP_CHEMKIN .EQ. 1 ) THEN
        DO i = 1, nx*ny*nz
! calculate heat capacity C_p and formation energy of mixture
           HEAT_CAPACITY    = C_0_R
           FORMATION_ENERGY = C_0_R
           WMEAN_INV        = C_0_R
           DO is = 1,NSP-1
              HEAT_CAPACITY    = HEAT_CAPACITY&
                   + s(i,is)*(THERMO_AI(1,1,is)-THERMO_AI(1,1,NSP))
              FORMATION_ENERGY = FORMATION_ENERGY&
                   + s(i,is)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
              WMEAN_INV        = WMEAN_INV&
                   + s(i,is)*(WGHT_INV(is)-WGHT_INV(NSP))
           ENDDO
           HEAT_CAPACITY    = HEAT_CAPACITY    + THERMO_AI(1,1,NSP)
           FORMATION_ENERGY = FORMATION_ENERGY + THERMO_AI(6,1,NSP)
           WMEAN_INV        = WMEAN_INV        + WGHT_INV(NSP)
! solve for T; go from C_p to C_v
           T(i)    = (e(i)-FORMATION_ENERGY)/(HEAT_CAPACITY-GRATIO*WMEAN_INV)
        ENDDO

! -------------------------------------------------------------------
! General case
! -------------------------------------------------------------------
     ELSE
        CALL TLAB_WRITE_ASCII(efile, 'THERMO_CALORIC_TEMPERATURE. General case undeveloped.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)

     ENDIF

  ENDIF

  RETURN
END SUBROUTINE THERMO_CALORIC_TEMPERATURE
