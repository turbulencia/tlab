#include "dns_const.h"

!########################################################################
!#
!# Calculate diagnostic variables
!#
!########################################################################
subroutine FI_DIAGNOSTIC(nx, ny, nz, q, s)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use TLAB_VARS, only: imode_eqns
    use TLAB_VARS, only: damkohler
    use TLab_Arrays, only: wrk1d, wrk3d
    use Thermodynamics
    use Gravity, only: buoyancy, bbackground
    use THERMO_THERMAL
    use THERMO_CALORIC
    use THERMO_AIRWATER
    use THERMO_ANELASTIC
    use Averages, only: AVG1V2D_V

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(INOUT) :: q(nx*ny*nz, inb_flow_array)
    real(wp), intent(INOUT) :: s(nx*ny*nz, inb_scal_array)

    ! -------------------------------------------------------------------
    ! ###################################################################
    select case (imode_eqns)
    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        if (buoyancy%type == EQNS_BOD_NORMALIZEDMEAN .or. buoyancy%type == EQNS_BOD_SUBTRACTMEAN) then      ! Calculate mean buoyancy background
            call AVG1V2D_V(nx, ny, nz, 1, s(:, 1), bbackground, wrk1d)
        end if
    end select

    select case (imode_thermo)
    case (THERMO_TYPE_ANELASTIC)
        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then       ! Calculate liquid content q_l
            call THERMO_ANELASTIC_PH(nx, ny, nz, s(:, 2), s(1, 1))

        end if

    case (THERMO_TYPE_LINEAR)
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then                             ! Calculate liquid content q_l
            call THERMO_AIRWATER_LINEAR(nx*ny*nz, s, s(:, inb_scal_array))

        end if

    case (THERMO_TYPE_COMPRESSIBLE)
#define e(j)    q(j,4)
#define rho(j)  q(j,5)
#define p(j)    q(j,6)
#define T(j)    q(j,7)
#define vis(j)  q(j,8)

        call THERMO_CALORIC_TEMPERATURE(nx*ny*nz, s, e(1), rho(1), T(1), wrk3d)
        call THERMO_THERMAL_PRESSURE(nx*ny*nz, s, rho(1), T(1), p(1))
        if (any([EQNS_TRANS_SUTHERLAND, EQNS_TRANS_POWERLAW] == itransport)) call THERMO_VISCOSITY(nx*ny*nz, T(1), vis(1))

    end select

    return
end subroutine FI_DIAGNOSTIC
