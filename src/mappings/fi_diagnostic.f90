#include "types.h"
#include "dns_const.h"

!########################################################################
!#
!# Calculate diagnostic variables
!#
!########################################################################
SUBROUTINE FI_DIAGNOSTIC( nx,ny,nz, q,s, wrk3d )

  USE TLAB_VARS,     ONLY : inb_flow_array, inb_scal_array
  USE TLAB_VARS,     ONLY : imode_eqns, itransport, damkohler
  USE TLAB_VARS,     ONLY : epbackground,pbackground
  USE THERMO_GLOBAL,  ONLY : imixture

  IMPLICIT NONE

  TINTEGER, INTENT(IN   ) :: nx,ny,nz
  TREAL,    INTENT(INOUT) :: q(nx*ny*nz,inb_flow_array)
  TREAL,    INTENT(INOUT) :: s(nx*ny*nz,inb_scal_array)
  TREAL,    INTENT(INOUT) :: wrk3d(nx*ny*nz)

  ! -------------------------------------------------------------------
  ! ###################################################################
  SELECT CASE( imode_eqns )
  CASE( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC )
    IF      ( imixture == MIXT_TYPE_AIRWATER .AND. damkohler(3) <= C_0_R ) THEN ! Calculate q_l
      CALL THERMO_AIRWATER_PH(nx,ny,nz, s(1,2), s(1,1), epbackground,pbackground)

    ELSE IF ( imixture == MIXT_TYPE_AIRWATER_LINEAR                      ) THEN
      CALL THERMO_AIRWATER_LINEAR(nx,ny,nz, s, s(1,inb_scal_array))

    ENDIF

  CASE( DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL )
#define e(j)    q(j,4)
#define rho(j)  q(j,5)
#define p(j)    q(j,6)
#define T(j)    q(j,7)
#define vis(j)  q(j,8)

    CALL THERMO_CALORIC_TEMPERATURE(nx,ny,nz, s, e(1), rho(1), T(1), wrk3d)
    CALL THERMO_THERMAL_PRESSURE(nx,ny,nz, s, rho(1), T(1), p(1))
    IF ( itransport == EQNS_TRANS_SUTHERLAND .OR. itransport == EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(nx,ny,nz, T(1), vis(1))

  END SELECT

  RETURN
END SUBROUTINE FI_DIAGNOSTIC
