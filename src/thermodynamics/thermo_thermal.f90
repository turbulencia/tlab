#include "dns_const.h"

! Implementation of thermal equation of state
! Default is NSP species represented by the first NSP-1 mass fractions in array s
! Airwater considers the special case where qt=qv+ql and ql instead of qv and ql are used as prognostic

! we use explicit-shape arrays in argements for the routines to be calleable with different array ranks

module THERMO_THERMAL
    use TLAB_CONSTANTS, only: wp, wi
    use THERMO_VARS, only: imixture, THERMO_R, MRATIO, RRATIO, NSP
    implicit none
    private

    integer(wi) ij, is
    real(wp) RMEAN

    public :: THERMO_THERMAL_DENSITY
    public :: THERMO_THERMAL_TEMPERATURE
    public :: THERMO_THERMAL_PRESSURE

contains
    !########################################################################
    ! Calculate rho from T, p and composition
    !########################################################################
    subroutine THERMO_THERMAL_DENSITY(ijmax, s, p, T, rho)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), p(ijmax), T(ijmax)
        real(wp), intent(out) :: rho(ijmax)

! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            rho(:) = MRATIO*p(:)/T(:)

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.
!     DO ij = 1,ijmax
!#define MACRO_ZINPUT s(ij,inb_scal)
!#include "dns_chem_mass.h"
!        rho(ij) = MRATIO*p(ij)*WMEAN/T(ij)
!     ENDDO

        case (MIXT_TYPE_AIRWATER)
            rho(:) = p(:)/(T(:)*(THERMO_R(2) + s(:, 1)*(THERMO_R(1) - THERMO_R(2)) - s(:, 2)*THERMO_R(1)))

        case default
            do ij = 1, ijmax
                RMEAN = THERMO_R(NSP)
                do is = 1, NSP - 1
                    RMEAN = RMEAN + s(ij, is)*(THERMO_R(is) - THERMO_R(NSP))
                end do
                rho(ij) = p(ij)/(RMEAN*T(ij))
            end do

        end select

    end subroutine THERMO_THERMAL_DENSITY

    !########################################################################
    ! Calculate T from rho, p and composition; same as before but swapping T and rho
    !########################################################################
    subroutine THERMO_THERMAL_TEMPERATURE(ijmax, s, p, rho, T)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), p(ijmax), rho(ijmax)
        real(wp), intent(out) :: T(ijmax)

! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            T(:) = MRATIO*p(:)/rho(:)

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.
!     DO ij = 1,ijmax
!#define MACRO_ZINPUT s(ij,inb_scal)
!#include "dns_chem_mass.h"
!        T(ij) = MRATIO*p(ij)*WMEAN/rho(ij)
!     ENDDO

        case (MIXT_TYPE_AIRWATER)
            T(:) = p(:)/(rho(:)*(THERMO_R(2) + s(:, 1)*(THERMO_R(1) - THERMO_R(2)) - s(:, 2)*THERMO_R(1)))

        case default
            do ij = 1, ijmax
                RMEAN = THERMO_R(NSP)
                do is = 1, NSP - 1
                    RMEAN = RMEAN + s(ij, is)*(THERMO_R(is) - THERMO_R(NSP))
                end do
                T(ij) = p(ij)/(RMEAN*rho(ij))
            end do

        end select

    end subroutine THERMO_THERMAL_TEMPERATURE

    !########################################################################
    ! Calculate T from rho, p and composition; same as before but swapping T and rho
    !########################################################################
    subroutine THERMO_THERMAL_PRESSURE(ijmax, s, rho, T, p)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), rho(ijmax), T(ijmax)
        real(wp), intent(out) :: p(ijmax)

! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            p(:) = RRATIO*rho(:)*T(:)

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.
!     DO ij = 1,ijmax
!#define MACRO_ZINPUT s(ij,inb_scal)
!#include "dns_chem_mass.h"
!        p(ij) = rho(ij)*T(ij)/(MRATIO*WMEAN)
!     ENDDO

        case (MIXT_TYPE_AIRWATER)
            p(:) = rho(:)*T(:)*(THERMO_R(2) + s(:, 1)*(THERMO_R(1) - THERMO_R(2)) - s(:, 2)*THERMO_R(1))

        case default
            do ij = 1, ijmax
                RMEAN = THERMO_R(NSP)
                do is = 1, NSP - 1
                    RMEAN = RMEAN + s(ij, is)*(THERMO_R(is) - THERMO_R(NSP))
                end do
                p(ij) = rho(ij)*T(ij)*RMEAN
            end do

        end select

    end subroutine THERMO_THERMAL_PRESSURE

end module THERMO_THERMAL
