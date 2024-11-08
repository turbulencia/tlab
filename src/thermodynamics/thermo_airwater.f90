#include "dns_const.h"

! Implementation of airwater cases
! Airwater considers the special case where qt=qv+ql and ql instead of qv and ql are used as prognostic

! we use explicit-shape arrays in arguments for the routines to be callable with different array ranks

module THERMO_AIRWATER
    use TLab_Constants, only: wp, wi, small_wp, big_wp
    use TLAB_VARS, only: inb_scal
    use Thermodynamics, only: imixture, CRATIO_INV, NCP, THERMO_AI, THERMO_TLIM
    use Thermodynamics, only: THERMO_PSAT, NPSAT, Thermo_Psat_Polynomial
    use Thermodynamics, only: Cd, Cdv, Cvl, Cdl, Ld, Lv, Ldv, Lvl, Ldl, Rd, Rdv, Rv, rd_ov_rv
    use Thermodynamics, only: dsmooth, NEWTONRAPHSON_ERROR
    use Thermodynamics, only: thermo_param
    use THERMO_THERMAL
    implicit none
    private

    integer(wi) ij, is, ipsat
    real(wp) RMEAN, CPMEAN, T_LOC, psat, qsat, dsmooth_loc, dummy
    integer(wi) inr, nrmax
    real(wp) B_LOC(10), FUN, DER, ERROR_LOC

    public :: THERMO_AIRWATER_PT
    public :: THERMO_AIRWATER_RP
    public :: THERMO_AIRWATER_PH_RE
    public :: THERMO_AIRWATER_RE
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

        real(wp) dqldqt

        ! ###################################################################
        call Thermo_Psat_Polynomial(ijmax, T, s(1, 2))
        do ij = 1, ijmax
            ! this is really the vapor content
            qsat = 1.0_wp/(p(ij)/s(ij, 2) - 1.0_wp)*rd_ov_rv*(1.0_wp - s(ij, 1))
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
    !# Calculate T and liquid content from rho, p and water content using thermal equation of state.
    !# The difference with THERMO_THERMAL_TEMPERATURE is that the equilibrium
    !# partition of q_t between vapor and liquid here is not known and q_l has to be calculated.
    !########################################################################
    subroutine THERMO_AIRWATER_RP(ijmax, s, p, rho, T, dqldqt)
        integer(wi) ijmax
        real(wp), intent(in) :: p(ijmax), rho(ijmax)
        real(wp), intent(inout) :: s(ijmax, *), dqldqt(ijmax)
        real(wp), intent(out) :: T(ijmax)

        ! -------------------------------------------------------------------
        real(wp) B_LOC_CONST_1, B_LOC_CONST_2, alpha

        ! ###################################################################
        ERROR_LOC = 0.0_wp
        ! maximum number of iterations in Newton-Raphson
        nrmax = 3

        ! reference case q_l = 0
        T(:) = p(:)/(rho(:)*(Rd + s(:, 1)*Rdv - s(:, 2)*Rv))

        ! -------------------------------------------------------------------
        ! calculate saturation specific humidity q_s(\rho,p)
        ! -------------------------------------------------------------------
        if (dsmooth <= 0.0_wp) then
            ! calculate saturation specific humidity, in array s(1,2).
            ! Thermo_Psat_Polynomial is duplicated here to avoid array calls
            do ij = 1, ijmax
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*T(ij) + THERMO_PSAT(ipsat)
                end do
                s(ij, 2) = psat/(rho(ij)*T(ij)*Rv)
            end do

        else
            ! initialize homogeneous data
            do ij = 1, 9
                B_LOC(ij) = THERMO_PSAT(ij)*(rd_ov_rv - 1.0_wp)
            end do
            B_LOC_CONST_1 = B_LOC(1)
            B_LOC_CONST_2 = B_LOC(2)

            ! loop on all points
            do ij = 1, ijmax
                B_LOC(1) = B_LOC_CONST_1 + p(ij)
                B_LOC(2) = B_LOC_CONST_2 - rho(ij)*Rd

                ! Newton-Raphson
                t_loc = T(ij)
                do inr = 1, nrmax
                    FUN = B_LOC(9)
                    DER = 0.0_wp
                    do is = 8, 1, -1
                        FUN = FUN*t_loc + B_LOC(is)
                        DER = DER*t_loc + B_LOC(is + 1)*real(is, wp)
                    end do
                    t_loc = t_loc - FUN/DER
                end do
                ERROR_LOC = max(ERROR_LOC, abs(FUN/DER)/t_loc)

                ! calculate saturation specific humidity, in array s(1,2).
                s(ij, 2) = (p(ij) - rho(ij)*Rd*t_loc)/(rho(ij)*Rdv*t_loc)

                ! calculate dqldqt
                qsat = s(ij, 2)
                alpha = -Lvl - (Cvl + CRATIO_INV*Rv)*t_loc
                alpha = alpha/(t_loc*CRATIO_INV*Rv) - 1.0_wp

                dummy = p(ij)/(qsat*rho(ij)*Rv*t_loc)

                dqldqt(ij) = 1.0_wp - alpha*rd_ov_rv/(1.0_wp + dummy)
            end do

        end if

        ! -------------------------------------------------------------------
        ! check if condesation occurs and, if so, solve nonlinear equation
        ! -------------------------------------------------------------------
        ! initialize homogeneous data
        B_LOC(1:9) = THERMO_PSAT(1:9)

        ! loop on all points
        do ij = 1, ijmax
            qsat = s(ij, 2)

            if (qsat > s(ij, 1)) then
                s(ij, 2) = 0.0_wp
                if (dsmooth > 0.0_wp) then
                    dsmooth_loc = dsmooth*qsat
                    s(ij, 2) = dsmooth_loc*dqldqt(ij) &
                               *log(exp((s(ij, 1) - qsat)/dsmooth_loc) + 1.0_wp)
                    ! change T consistently
                    T(ij) = p(ij)/(((s(ij, 1) - s(ij, 2))*Rv + (1.0_wp - s(ij, 1))*Rd)*rho(ij))
                end if

                ! if q_s < q_t, then we have to repeat calculation of T
            else
                B_LOC(1) = THERMO_PSAT(1) - p(ij)
                B_LOC(2) = THERMO_PSAT(2) + (1.0_wp - s(ij, 1))*rho(ij)*Rd

                ! Newton-Raphson
                do inr = 1, nrmax
                    FUN = B_LOC(9)
                    DER = 0.0_wp
                    do is = 8, 1, -1
                        FUN = FUN*T(ij) + B_LOC(is)
                        DER = DER*T(ij) + B_LOC(is + 1)*real(is, wp)
                    end do
                    T(ij) = T(ij) - FUN/DER
                end do
                ERROR_LOC = max(ERROR_LOC, abs(FUN/DER)/T(ij))

                ! calculate saturation specific humidity, in array s(1,2).
                ! Thermo_Psat_Polynomial is duplicated here to avoid array calls
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*T(ij) + THERMO_PSAT(ipsat)
                end do
                s(ij, 2) = psat/(rho(ij)*T(ij)*Rv)

                ! liquid content
                s(ij, 2) = s(ij, 1) - s(ij, 2)
                if (dsmooth > 0.0_wp) then
                    dsmooth_loc = dsmooth*qsat
                    s(ij, 2) = dsmooth_loc*dqldqt(ij) &
                               *log(exp((s(ij, 1) - qsat)/dsmooth_loc) + 1.0_wp) &
                               + s(ij, 2) - dqldqt(ij)*(s(ij, 1) - qsat)
                    ! change T consistently
                    T(ij) = p(ij)/ &
                            (((s(ij, 1) - s(ij, 2))*Rv + (1.0_wp - s(ij, 1))*Rd)*rho(ij))
                end if
            end if

        end do

        return
    end subroutine THERMO_AIRWATER_RP

    !########################################################################
    !Calculating the equilibrium T and q_l for given enthalpy and pressure
    !Iterative method based on (rho,e,q_i)
    !########################################################################
    subroutine THERMO_AIRWATER_PH_RE(ijmax, s, p, h, T)
        integer(wi) ijmax
        real(wp), intent(in) :: p(ijmax), h(ijmax)
        real(wp), intent(inout) :: s(ijmax, *)
        real(wp), intent(out) :: T(ijmax)

        ! -------------------------------------------------------------------
        integer(wi) iter, niter
        integer(wi) ij_loc
        real(wp) r_loc(1), e_loc(1), t_loc(1), s_loc(2), dummy(1)

        ! ###################################################################
        niter = 5

        do ij_loc = 1, ijmax
            ! initialize, q_l=0
            s_loc(1) = s(ij_loc, 1)
            s_loc(2) = 0.0_wp
            t_loc(1) = (h(ij_loc) - Ld - s(ij_loc, 1)*Ldv)/(Cd + s(ij_loc, 1)*Cdv)

            do iter = 1, niter      ! iteration
                ! calculate density from temperature/composition
                call THERMO_THERMAL_DENSITY(1, s_loc, p(ij_loc), t_loc, r_loc)
                ! calculate energy
                e_loc = h(ij_loc) - CRATIO_INV*p(ij_loc)/r_loc(1)
                ! solve equilibrium (rho,e,q_i)
                call THERMO_AIRWATER_RE(1, s_loc, e_loc, r_loc, t_loc(:), dummy(:))
            end do
            s(ij_loc, 2) = s_loc(2)
            T(ij_loc) = t_loc(1)

        end do

        return
    end subroutine THERMO_AIRWATER_PH_RE

!########################################################################
!# Calculate T and liquid content from rho, e and water content using caloric equation of state.
!########################################################################
    subroutine THERMO_AIRWATER_RE(ijmax, s, e, rho, T, dqldqt)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: e(ijmax), rho(ijmax)
        real(wp), intent(out) :: T(ijmax)
        real(wp), intent(inout) :: s(ijmax, *), dqldqt(ijmax)

        ! -------------------------------------------------------------------
        integer(wi) i
        real(wp) HEAT_CAPACITY_LD, HEAT_CAPACITY_LV, HEAT_CAPACITY_VD
        real(wp) B_LOC_CONST_2, B_LOC_CONST_3
        real(wp) alpha, dummy1, dummy2

        ! ###################################################################
        NEWTONRAPHSON_ERROR = 0.0_wp
        ! maximum number of iterations in Newton-Raphson
        nrmax = 3

        ! reference case q_l = 0
        do ij = 1, ijmax
            T(ij) = e(ij) - (Ld + s(ij, 1)*Ldv)
            CPMEAN = Cd + s(ij, 1)*Cdv
            RMEAN = Rd + s(ij, 1)*Rdv
            T(ij) = (e(ij) - (Ld + s(ij, 1)*Ldv))/(CPMEAN - RMEAN*CRATIO_INV)
        end do

        ! -------------------------------------------------------------------
        ! calculate saturation specific humidity q_s(\rho,e)
        ! -------------------------------------------------------------------
        if (dsmooth <= 0.0_wp) then
            ! calculate saturation specific humidity, in array s(1,2).
            ! Thermo_Psat_Polynomial is duplicated here to avoid array calls
            do i = 1, ijmax
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*T(i) + THERMO_PSAT(ipsat)
                end do
                s(i, 2) = psat/(rho(i)*T(i)*Rv)
            end do

        else
            ! initialize homogeneous data
            HEAT_CAPACITY_VD = Rdv*CRATIO_INV - Cdv
            do i = 1, 9
                B_LOC(i) = -THERMO_PSAT(i)*Ldv
            end do
            B_LOC(10) = 0.0_wp
            do i = 2, 10
                B_LOC(i) = B_LOC(i) + THERMO_PSAT(i - 1)*HEAT_CAPACITY_VD
            end do
            B_LOC_CONST_2 = B_LOC(2)
            B_LOC_CONST_3 = B_LOC(3)

            ! loop on all points
            do i = 1, ijmax
                B_LOC(2) = B_LOC_CONST_2 + rho(i)*Rv*(e(i) - Lv)
                B_LOC(3) = B_LOC_CONST_3 - rho(i)*Rv*(Cd - CRATIO_INV*Rd)
                ! Newton-Raphson
                t_loc = T(i)
                do inr = 1, nrmax
                    FUN = B_LOC(10)
                    DER = 0.0_wp
                    do is = 9, 1, -1
                        FUN = FUN*t_loc + B_LOC(is)
                        DER = DER*t_loc + B_LOC(is + 1)*real(is, wp)
                    end do
                    t_loc = t_loc - FUN/DER
                end do
                NEWTONRAPHSON_ERROR = max(NEWTONRAPHSON_ERROR, abs(FUN/DER)/t_loc)
                ! calculate saturation specific humidity, in array s(1,2).
                ! Thermo_Psat_Polynomial is duplicated here to avoid routine calls
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*t_loc + THERMO_PSAT(ipsat)
                end do
                s(i, 2) = psat/(rho(i)*t_loc*Rv)

                ! calculate dqldqt
                qsat = s(i, 2)
                alpha = -Lvl - (Cvl + CRATIO_INV*Rv)*t_loc
                alpha = alpha/(t_loc*CRATIO_INV*Rv) - 1.0_wp

                dummy1 = -Ldl - (Cdl + CRATIO_INV*Rd)*t_loc
                dummy1 = dummy1*qsat*alpha

                dummy2 = -Ldv + t_loc*CRATIO_INV*Rv*alpha*(alpha + 1.0_wp)
                dummy2 = dummy2*qsat + e(i) - Lv

                dqldqt(i) = 1.0_wp - dummy1/dummy2
            end do

        end if

        ! -------------------------------------------------------------------
        ! calculate final T and q_l
        ! -------------------------------------------------------------------
        ! initialize homogeneous data
        HEAT_CAPACITY_LV = Cvl + CRATIO_INV*Rv
        HEAT_CAPACITY_LD = Cdl + CRATIO_INV*Rd
        HEAT_CAPACITY_VD = HEAT_CAPACITY_LD - HEAT_CAPACITY_LV
        do i = 1, 9
            B_LOC(i) = THERMO_PSAT(i)*Lvl
        end do
        B_LOC(10) = 0.0_wp
        do i = 2, 10
            B_LOC(i) = B_LOC(i) + THERMO_PSAT(i - 1)*HEAT_CAPACITY_LV
        end do
        B_LOC_CONST_2 = B_LOC(2)
        B_LOC_CONST_3 = B_LOC(3)

        ! loop on all points
        do i = 1, ijmax
            qsat = s(i, 2)

            if (qsat >= s(i, 1)) then
                s(i, 2) = 0.0_wp
                if (dsmooth > 0.0_wp) then
                    dsmooth_loc = dsmooth*qsat
                    s(i, 2) = dsmooth_loc*dqldqt(i) &
                              *log(exp((s(i, 1) - qsat)/dsmooth_loc) + 1.0_wp)
                    ! change T consistently
                    T(i) = (e(i) - s(i, 1)*Ldv - Lv - s(i, 2)*Lvl)/ &
                           (s(i, 1)*HEAT_CAPACITY_VD + Cd - CRATIO_INV*Rd + s(i, 2)*HEAT_CAPACITY_LV)
                end if

                ! if q_s < q_t, then we have to repeat calculation of T
            else
                B_LOC(2) = B_LOC_CONST_2 + rho(i)*Rv*(e(i) - s(i, 1)*Ldl - Ld)
                B_LOC(3) = B_LOC_CONST_3 - rho(i)*Rv*(s(i, 1)*HEAT_CAPACITY_LD + Cd - CRATIO_INV*Rd)
                ! IF ( dsmooth .GT. 0.0_wp ) THEN
                ! dsmooth_loc = dsmooth*qsat
                ! alpha =-( dsmooth_loc*LOG(EXP((s(i,1)-qsat)/dsmooth_loc)+1.0_wp)
                ! $                 -(s(i,1)-qsat) )*dqldqt(i)
                ! B_LOC(2) = B_LOC(2) + rho(i)*Rv*alpha*Lvl
                ! B_LOC(3) = B_LOC(3) + rho(i)*Rv*alpha*HEAT_CAPACITY_LV
                ! ENDIF

                ! Newton-Raphson
                do inr = 1, nrmax
                    FUN = B_LOC(10)
                    DER = 0.0_wp
                    do is = 9, 1, -1
                        FUN = FUN*T(i) + B_LOC(is)
                        DER = DER*T(i) + B_LOC(is + 1)*real(is, wp)
                    end do
                    T(i) = T(i) - FUN/DER
                end do
                NEWTONRAPHSON_ERROR = max(NEWTONRAPHSON_ERROR, abs(FUN/DER)/T(i))

                ! calculate saturation specific humidity, in array s(1,2).
                ! Thermo_Psat_Polynomial is duplicated here to avoid routine calls
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*T(i) + THERMO_PSAT(ipsat)
                end do
                s(i, 2) = psat/(rho(i)*T(i)*Rv)

                ! liquid content
                s(i, 2) = s(i, 1) - s(i, 2)
                if (dsmooth > 0.0_wp) then
                    dsmooth_loc = dsmooth*qsat
                    s(i, 2) = dsmooth_loc*dqldqt(i) &
                              *log(exp((s(i, 1) - qsat)/dsmooth_loc) + 1.0_wp) &
                              + s(i, 2) - dqldqt(i)*(s(i, 1) - qsat)
                    ! change T consistently
                    T(i) = (e(i) - s(i, 1)*Ldv - Lv - s(i, 2)*Lvl)/ &
                           (s(i, 1)*HEAT_CAPACITY_VD + Cd - CRATIO_INV*Rd &
                            + s(i, 2)*HEAT_CAPACITY_LV)
                end if
            end if

        end do

        return
    end subroutine THERMO_AIRWATER_RE

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
