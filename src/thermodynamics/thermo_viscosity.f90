#include "dns_const.h"

!########################################################################
!# Calculate the nondimensional dynamic viscosity \mu=\mu(T) depending
!# on the given flag itransport:
!########################################################################
subroutine THERMO_VISCOSITY(ijmax, T, mu)
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: itransport
    implicit none

    integer(wi) ijmax
    real(wp), intent(in):: T(ijmax)
    real(wp), intent(out):: mu(ijmax)

! #######################################################################
    select case (itransport)

    case (EQNS_TRANS_SUTHERLAND)    ! not yet implemented
        mu(:) = 1.0_wp

    case (EQNS_TRANS_POWERLAW)
        mu(:) = T(:)**0.7_wp

    case default
        mu(:) = 1.0_wp
    end select

    return
end subroutine THERMO_VISCOSITY
