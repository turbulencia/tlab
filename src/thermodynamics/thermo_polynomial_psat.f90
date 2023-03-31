#include "dns_error.h"

subroutine THERMO_POLYNOMIAL_PSAT(ijmax, T, p)
    use TLAB_CONSTANTS, only: wp, wi
    use THERMO_VARS, only: THERMO_PSAT, NPSAT

    implicit none

    integer(wi) ijmax
    real(wp) T(ijmax)
    real(wp) p(ijmax)

! -------------------------------------------------------------------
    integer(wi) i, ipsat

! ###################################################################
    if (NPSAT > 0) then
        do i = 1, ijmax
            p(i) = THERMO_PSAT(NPSAT)
            do ipsat = NPSAT - 1, 1, -1
                p(i) = p(i)*T(i) + THERMO_PSAT(ipsat)
            end do
        end do
    else
        p(:) = 0.0_wp
    end if

    return
end subroutine THERMO_POLYNOMIAL_PSAT

! ###################################################################
! ###################################################################
subroutine THERMO_POLYNOMIAL_DPSAT(ijmax, T, dp)
    use TLAB_CONSTANTS, only: wp, wi
    use THERMO_VARS, only: THERMO_PSAT, NPSAT

    implicit none

    integer(wi) ijmax
    real(wp) T(ijmax)
    real(wp) dp(ijmax)

! -------------------------------------------------------------------
    integer(wi) i, ipsat

! ###################################################################
    if (NPSAT > 0) then
        do i = 1, ijmax
            dp(i) = 0.0_wp
            do ipsat = NPSAT - 1, 1, -1
                dp(i) = dp(i)*T(i) + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
            end do
        end do
    else
        dp(:) = 0.0_wp
    end if

    return
end subroutine THERMO_POLYNOMIAL_DPSAT
