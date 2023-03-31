#include "dns_const.h"

! Implementation of airwater cases
! Airwater considers the special case where qt=qv+ql and ql instead of qv and ql are used as prognostic

! we use explicit-shape arrays in arguments for the routines to be callable with different array ranks

module THERMO_AIRWATER
    use TLAB_VARS, only: inb_scal
    use TLAB_CONSTANTS, only: wp, wi, small_wp, big_wp
    use THERMO_VARS, only: imixture, MRATIO, RRATIO, GRATIO, CRATIO_INV, THERMO_PSAT, NPSAT
    use THERMO_VARS, only: imixture, Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl, Cdl, Cl, rd_ov_rv, rd_ov_cd, PREF_THETA
    use THERMO_VARS, only: dsmooth
    use THERMO_VARS, only: thermo_param
    implicit none
    private

    integer(wi) ij, is
    real(wp) RMEAN, qsat, dqldqt, dsmooth_loc, dummy

    public :: THERMO_AIRWATER_PT
    public :: THERMO_AIRWATER_LINEAR
    public :: THERMO_AIRWATER_LINEAR_SOURCE

contains
    !########################################################################
    !# Calculate liquid content from p, T and water content.
    !########################################################################
    subroutine THERMO_AIRWATER_PT(ijmax, s, p, T)
        integer(wi) ijmax
        real(wp), intent(in) :: p(ijmax), T(ijmax)
        real(wp), intent(out) :: s(ijmax, *)

        ! ###################################################################
        call THERMO_POLYNOMIAL_PSAT(ijmax, T, s(1, 2))
        do ij = 1, ijmax
            ! this is really the vapor content
            qsat = 1.0_wp/(MRATIO*p(ij)/s(ij, 2) - 1.0_wp)*rd_ov_rv*(1 - s(ij, 1))
            if (qsat >= s(ij, 1)) then
                s(ij, 2) = 0.0_wp
            else
                s(ij, 2) = s(ij, 1) - qsat
            end if
            if (dsmooth > 0.0_wp) then
                qsat = qsat/(1 - s(ij, 1))
                dqldqt = 1.0_wp + qsat
                ! this is the real qsat
                qsat = qsat/(1.0_wp + qsat)
                dsmooth_loc = dsmooth*qsat
                s(ij, 2) = dsmooth_loc*dqldqt &
                           *log(exp((s(ij, 1) - qsat)/dsmooth_loc) + 1.0_wp)
            end if
        end do

        return
    end subroutine THERMO_AIRWATER_PT

    !########################################################################
    !# Calculating the equilibrium liquid from the mixture fraction and the
    !# enthalpy deviation according to the linearized equilibrium thermodynamics
    !########################################################################
    subroutine THERMO_AIRWATER_LINEAR(ijmax, s, l)
        integer(wi), intent(IN) :: ijmax
        real(wp), dimension(ijmax, inb_scal), intent(IN) :: s     ! chi, psi
        real(wp), dimension(ijmax), intent(OUT) :: l     ! normalized liquid

        ! -------------------------------------------------------------------
        real(wp) dummy2

        ! ###################################################################
        ! Calculating \xi
        if (inb_scal == 1) then
            l = 1.0_wp + thermo_param(1)*s(:, 1)
        else
            l = 1.0_wp + thermo_param(1)*s(:, 1) + thermo_param(2)*s(:, 2)
        end if

        ! Calculating liquid; \xi is overwritten in this routine
        if (abs(thermo_param(inb_scal + 1)) < small_wp) then
            do ij = 1, ijmax
                l(ij) = max(l(ij), 0.0_wp)
            end do

        else
            dummy = thermo_param(inb_scal + 1)
            dummy2 = 1.0_wp/dummy
            !     l = dummy *LOG( EXP(dummy2 *l) +1.0_wp )
            do ij = 1, ijmax
                l(ij) = dummy*log(exp(dummy2*l(ij)) + 1.0_wp)
            end do

        end if

        return
    end subroutine THERMO_AIRWATER_LINEAR

    !########################################################################
    !########################################################################
    subroutine THERMO_AIRWATER_LINEAR_SOURCE(ijmax, s, xi, der1, der2)
        integer(wi), intent(IN) :: ijmax
        real(wp), intent(IN) :: s(ijmax, inb_scal)        ! chi, psi
        real(wp), intent(OUT) :: xi(ijmax), der1(ijmax), der2(ijmax)

        ! ###################################################################
        if (inb_scal == 1) then
            xi(:) = 1.0_wp + thermo_param(1)*s(:, 1)
        else
            xi(:) = 1.0_wp + thermo_param(1)*s(:, 1) + thermo_param(2)*s(:, 2)
        end if

        if (abs(thermo_param(inb_scal + 1)) < small_wp) then
            der2(:) = big_wp

            do ij = 1, ijmax
                if (xi(ij) <= 0.0_wp) then
                    der1(ij) = 0.0_wp
                else
                    der1(ij) = 1.0_wp
                end if
            end do

        else
            ! Formulation in terms of the exponential might lead to NAN because of too large numbers!
            !     dummy=-1.0_wp/thermo_param(inb_scal+1)
            !     der1 = 1.0_wp / ( 1.0_wp + EXP(dummy *xi) )

            !     der2 = (der1-1.0_wp) *der1 *dummy

            dummy = 0.5_wp/thermo_param(inb_scal + 1)
            der1(:) = 0.5_wp*(1.0_wp + tanh(dummy*xi(:)))

            dummy = 1.0_wp/thermo_param(inb_scal + 1)
            der2(:) = (1.0_wp - der1(:))*der1(:)*dummy

        end if

        return
    end subroutine THERMO_AIRWATER_LINEAR_SOURCE

end module THERMO_AIRWATER
