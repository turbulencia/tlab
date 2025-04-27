#include "dns_const.h"

!########################################################################
!#
!# Calculate thermodynamic properties of airwater in the anelastic approximation,
!# when a background state is given in the form of appropriate profiles
!#
!# s1 is specific static energy
!# s2 is total water specific humidity
!# s3 is liquid water specific humidity, if any
!#
!########################################################################

module Thermo_Anelastic
    use TLab_Constants, only: wp, wi
    use Thermodynamics, only: imixture
    use Thermodynamics, only: GRATIO, scaleheightinv   ! anelastic parameters
    use Thermodynamics, only: THERMO_PSAT, NPSAT, Thermo_Psat_Polynomial
    use Thermodynamics, only: Rv, Rd, Rdv, Cd, Cdv, Lv0, Ld, Ldv, Cvl, Cdl, Cl, rd_ov_rv, PREF_1000
    use Thermodynamics, only: MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER
    implicit none
    private

    integer(wi) ij, i, j, jk, is, ipsat
    real(wp) RMEAN, P_LOC, E_LOC, T_LOC, R_LOC, R_LOC_INV, psat

    ! real(wp), public :: scaleheight, GRATIO

    public :: Thermo_Anelastic_PH
    public :: Thermo_Anelastic_TEMPERATURE
    public :: Thermo_Anelastic_STATIC_L
    public :: Thermo_Anelastic_BUOYANCY
    public :: Thermo_Anelastic_WEIGHT_ADD
    public :: Thermo_Anelastic_WEIGHT_INPLACE
    public :: Thermo_Anelastic_WEIGHT_OUTPLACE
    public :: Thermo_Anelastic_WEIGHT_SUBTRACT

    public :: Thermo_Anelastic_DENSITY
    public :: Thermo_Anelastic_THETA
    public :: Thermo_Anelastic_THETA_V
    public :: Thermo_Anelastic_THETA_E
    public :: Thermo_Anelastic_THETA_L
    public :: Thermo_Anelastic_LAPSE_EQU
    public :: Thermo_Anelastic_LAPSE_FR
    public :: Thermo_Anelastic_VAPOR_PRESSURE
    public :: Thermo_Anelastic_DEWPOINT
    public :: Thermo_Anelastic_RELATIVEHUMIDITY
    public :: Thermo_Anelastic_STATIC_CONSTANTCP
    ! public :: Thermo_Anelastic_PH_RE

    ! background, reference profiles
    real(wp), allocatable, public :: epbackground(:)                    ! Potential energy
    real(wp), allocatable, public :: pbackground(:)                     ! Pressure background profile
    real(wp), allocatable, public :: tbackground(:)                     ! Temperature
    real(wp), allocatable, public :: rbackground(:), ribackground(:)    ! Density and its inverse

contains

    !########################################################################
    !########################################################################
    ! Kernel routines. They need to be fast.

    !########################################################################
    !# Calculating the equilibrium T and q_l for given enthalpy and pressure.
    !# Assumes often that THERMO_AI(6,1,1) = THERMO_AI(6,1,2) = 0
    !#
    !# Routine Thermo_Psat_Polynomial is duplicated here to avoid array calls
    !#
    !# Smoothing according to Eq. 25 in Mellado et al., TCFD, 2010
    !#
    !# s1 is total water specific humidity
    !# s2 is liquid water specific humidity
    !#
    !########################################################################
    subroutine Thermo_Anelastic_PH(nx, ny, nz, s, h)
        use Thermodynamics, only: dsmooth, NEWTONRAPHSON_ERROR

        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(inout) :: s(nx*ny*nz, *)
        real(wp), intent(in) :: h(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) inr, nrmax
        real(wp) ALPHA_1, ALPHA_2, BETA_1, BETA_2, alpha, beta, C_LN2_L
        real(wp) qsat, H_LOC, B_LOC(10), FUN, DER
        real(wp) dsmooth_loc, dqldqt, dqsdt!, qvequ
        real(wp) dummy

        ! ###################################################################
        NEWTONRAPHSON_ERROR = 0.0_wp
        ! maximum number of iterations in Newton-Raphson
        nrmax = 5

        ALPHA_1 = rd_ov_rv*Lv0
        ALPHA_2 = Lv0*(1.0_wp - rd_ov_rv)
        BETA_1 = rd_ov_rv*Cvl + Cd
        BETA_2 = Cdl - rd_ov_rv*Cvl

        C_LN2_L = log(2.0_wp)

        ! ###################################################################
        ij = 0
        do jk = 0, ny*nz - 1
            is = mod(jk, ny) + 1
            P_LOC = pbackground(is)
            E_LOC = epbackground(is)

            do i = 1, nx
                ij = ij + 1

                H_LOC = h(ij) - E_LOC ! enthalpy

                ! -------------------------------------------------------------------
                ! reference case assuming ql = 0
                ! -------------------------------------------------------------------
                s(ij, 2) = 0.0_wp
                T_LOC = H_LOC/(Cd + s(ij, 1)*Cdv)

                ! calculate saturation specific humidity q_sat(T,p)
                psat = THERMO_PSAT(NPSAT)
                do ipsat = NPSAT - 1, 1, -1
                    psat = psat*T_LOC + THERMO_PSAT(ipsat)
                end do
                dummy = rd_ov_rv/(P_LOC/psat - 1.0_wp)
                qsat = dummy/(1.0_wp + dummy)

                ! -------------------------------------------------------------------
                ! calculate smoothed piecewise linear contribution, if needed
                ! -------------------------------------------------------------------
                if (dsmooth > 0.0_wp) then
                    ! calculate dqsdt from dpsatdt (qs here mean qvequ)
                    dqsdt = 0.0_wp
                    do ipsat = NPSAT - 1, 1, -1
                        dqsdt = dqsdt*T_LOC + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                    end do
                    dqsdt = qsat/psat/(1.0_wp - psat/P_LOC)*dqsdt
                    ! calculate dqldqt
                    dqsdt = dqsdt/(Cd + qsat*Cdv)
                    dqldqt = (1.0_wp/(1.0_wp - qsat) + Cdv*T_LOC*dqsdt)/ &
                             (1.0_wp + (Lv0 - Cvl*T_LOC)*dqsdt)

                    dsmooth_loc = dsmooth*qsat
                    if (s(ij, 1) - qsat < 0.0_wp) then
                        s(ij, 2) = dqldqt*dsmooth_loc*log(exp((s(ij, 1) - qsat)/dsmooth_loc) + 1.0_wp)
                    else
                        s(ij, 2) = dqldqt*((s(ij, 1) - qsat) &
                                           + dsmooth_loc*(C_LN2_L - log(tanh((s(ij, 1) - qsat)/(2.0_wp*dsmooth_loc)) + 1.0_wp)))
                    end if

                end if

                ! -------------------------------------------------------------------
                ! if q_s < q_t, then we have to recalculate T
                ! -------------------------------------------------------------------
                if (qsat < s(ij, 1)) then
                    ! preparing Newton-Raphson
                    alpha = (ALPHA_1 + s(ij, 1)*ALPHA_2 + H_LOC)/P_LOC
                    beta = (BETA_1 + s(ij, 1)*BETA_2)/P_LOC
                    B_LOC(1) = H_LOC + s(ij, 1)*Lv0 - THERMO_PSAT(1)*alpha
                    do is = 2, 9
                        B_LOC(is) = THERMO_PSAT(is - 1)*beta - THERMO_PSAT(is)*alpha
                    end do
                    B_LOC(2) = B_LOC(2) - Cd - s(ij, 1)*Cdl
                    B_LOC(10) = THERMO_PSAT(9)*beta

                    ! executing Newton-Raphson
                    do inr = 1, nrmax
                        FUN = B_LOC(10)
                        DER = 0.0_wp
                        do is = 9, 1, -1
                            FUN = FUN*T_LOC + B_LOC(is)
                            DER = DER*T_LOC + B_LOC(is + 1)*real(is, wp)
                        end do
                        T_LOC = T_LOC - FUN/DER
                    end do
                    NEWTONRAPHSON_ERROR = max(NEWTONRAPHSON_ERROR, abs(FUN/DER)/T_LOC)

                    ! calculate equilibrium vapor specific humidity
                    psat = 0.0_wp
                    do ipsat = NPSAT, 1, -1
                        psat = psat*T_LOC + THERMO_PSAT(ipsat)
                    end do
                    !           dummy = rd_ov_rv /( P_LOC/psat -1.0_wp )
                    !           qvequ = dummy *( 1.0_wp -s(ij,1) )

                    if (dsmooth > 0.0_wp) then ! add correction
                        !              s(ij,2) = s(ij,2) +s(ij,1) -qvequ - (s(ij,1) -qsat) *dqldqt
                        s(ij, 2) = s(ij, 2) + s(ij, 1) - rd_ov_rv/(P_LOC/psat - 1.0_wp)*(1.0_wp - s(ij, 1)) - (s(ij, 1) - qsat)*dqldqt
                    else                           ! or calculate new
                        !              s(ij,2) =          s(ij,1) -qvequ
                        s(ij, 2) = s(ij, 1) - rd_ov_rv/(P_LOC/psat - 1.0_wp)*(1.0_wp - s(ij, 1))
                    end if

                end if

            end do
        end do

        return
    end subroutine Thermo_Anelastic_PH

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: T(nx*ny*nz)

! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_AIR)
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T(ij) = s(ij, 1) - E_LOC
                end do

            end do

        case (MIXT_TYPE_AIRVAPOR)
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T(ij) = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                end do

            end do

        case (MIXT_TYPE_AIRWATER)
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T(ij) = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                end do

            end do

        end select

        return
    end subroutine Thermo_Anelastic_TEMPERATURE

!########################################################################
! Calculating h_l - h; very similar to the temperature routine
!########################################################################
    subroutine Thermo_Anelastic_STATIC_L(nx, ny, nz, s, result)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: result(nx*ny*nz)

! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    result(ij) = s(ij, 1) - E_LOC
                    result(ij) = Cl*result(ij) + E_LOC - Lv0 - s(ij, 1)
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    result(ij) = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    result(ij) = Cl*result(ij) + E_LOC - Lv0 - s(ij, 1)
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    result(ij) = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    result(ij) = Cl*result(ij) + E_LOC - Lv0 - s(ij, 1)
                end do

            end do
        end if

        return
    end subroutine Thermo_Anelastic_STATIC_L

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_BUOYANCY(nx, ny, nz, s, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: b(nx*ny*nz)

! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)
                R_LOC = rbackground(is)
                R_LOC_INV = 1.0_wp/R_LOC

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    b(ij) = R_LOC_INV*(R_LOC - P_LOC/T_LOC)
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)
                R_LOC = rbackground(is)
                R_LOC_INV = 1.0_wp/R_LOC

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    RMEAN = Rd + s(ij, 2)*Rdv
                    b(ij) = R_LOC_INV*(R_LOC - P_LOC/(RMEAN*T_LOC))
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)
                R_LOC = rbackground(is)
                R_LOC_INV = 1.0_wp/R_LOC

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    RMEAN = Rd + s(ij, 2)*Rdv - s(ij, 3)*Rv
                    b(ij) = R_LOC_INV*(R_LOC - P_LOC/(RMEAN*T_LOC))
                end do

            end do
        end if

        return
    end subroutine Thermo_Anelastic_BUOYANCY

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_WEIGHT_INPLACE(nx, ny, nz, weight, a)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(*)
        real(wp), intent(inout) :: a(nx, ny*nz)

! ###################################################################
        do jk = 1, ny*nz
            j = mod(jk - 1, ny) + 1

            a(1:nx, jk) = a(1:nx, jk)*weight(j)

        end do

        return
    end subroutine Thermo_Anelastic_WEIGHT_INPLACE

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_WEIGHT_OUTPLACE(nx, ny, nz, weight, a, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(*)
        real(wp), intent(in) :: a(nx, ny*nz)
        real(wp), intent(out) :: b(nx, ny*nz)

! ###################################################################
        do jk = 1, ny*nz
            j = mod(jk - 1, ny) + 1

            b(1:nx, jk) = a(1:nx, jk)*weight(j)

        end do

        return
    end subroutine Thermo_Anelastic_WEIGHT_OUTPLACE

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_WEIGHT_ADD(nx, ny, nz, weight, a, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(*)
        real(wp), intent(in) :: a(nx, ny*nz)
        real(wp), intent(inout) :: b(nx, ny*nz)

! ###################################################################
        do jk = 1, ny*nz
            j = mod(jk - 1, ny) + 1

            b(1:nx, jk) = b(1:nx, jk) + a(1:nx, jk)*weight(j)

        end do

        return
    end subroutine Thermo_Anelastic_WEIGHT_ADD

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_WEIGHT_SUBTRACT(nx, ny, nz, weight, a, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(*)
        real(wp), intent(in) :: a(nx, ny*nz)
        real(wp), intent(inout) :: b(nx, ny*nz)

! ###################################################################
        do jk = 1, ny*nz
            j = mod(jk - 1, ny) + 1

            b(1:nx, jk) = b(1:nx, jk) - a(1:nx, jk)*weight(j)

        end do

        return
    end subroutine Thermo_Anelastic_WEIGHT_SUBTRACT

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_PRESSURE(nx, ny, nz, p)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(out) :: p(nx*ny*nz)

        ! ###################################################################
        ij = 0
        do jk = 0, ny*nz - 1
            is = mod(jk, ny) + 1

            do i = 1, nx
                ij = ij + 1
                p(ij) = pbackground(is)
            end do

        end do

        return
    end subroutine Thermo_Anelastic_PRESSURE

    !########################################################################
    !########################################################################
    ! Procedures for post-processing and initialization.
    ! We could write them explicitly as above using loops and only 1 array.
    ! But we favor clarity over speed.
    ! We write them in terms of a few routines despite requiring more memory and time

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_DENSITY(nx, ny, nz, s, rho, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: rho(nx*ny*nz), T(nx*ny*nz)

        ! ###################################################################
#define p rho

        call Thermo_Anelastic_PRESSURE(nx, ny, nz, p)
        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            rho = p/(T*Rd)

        case (MIXT_TYPE_AIRVAPOR)
            rho = p/(T*(Rd + s(:, 2)*Rdv))

        case (MIXT_TYPE_AIRWATER)
            rho = p/(T*(Rd + s(:, 2)*Rdv - s(:, 3)*Rv))

        end select

#undef p

        return
    end subroutine Thermo_Anelastic_DENSITY

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, pi)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(out) :: pi(nx*ny*nz)

        real(wp) kappa

        ! ###################################################################
        kappa = Rd/Cd*GRATIO

#define p pi

        call Thermo_Anelastic_PRESSURE(nx, ny, nz, p)

        pi = (PREF_1000/p)**kappa

#undef p

        return
    end subroutine Thermo_Anelastic_ONE_OV_EXNER

    !########################################################################
    ! Dry potential temperature
    !########################################################################
    subroutine Thermo_Anelastic_THETA(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        ! ###################################################################
#define locPi theta

        call Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, locPi)
        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)

        theta = T*locPi

#undef locPi

        return
    end subroutine Thermo_Anelastic_THETA

    !########################################################################
    ! Virtual Potential temperature
    !########################################################################
    subroutine Thermo_Anelastic_THETA_V(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        real(wp) Rdv_ov_Rd, Rv_ov_Rd

        ! ###################################################################
        Rdv_ov_Rd = Rdv/Rd
        Rv_ov_Rd = Rv/Rd

        call Thermo_Anelastic_THETA(nx, ny, nz, s, theta, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)

        case (MIXT_TYPE_AIRVAPOR)
            theta = theta*(1.0_wp + s(:, 2)*Rdv_ov_Rd)

        case (MIXT_TYPE_AIRWATER)
            theta = theta*(1.0_wp + s(:, 2)*Rdv_ov_Rd - s(:, 3)*Rv_ov_Rd)

        end select

        return
    end subroutine Thermo_Anelastic_THETA_V

    !########################################################################
    ! Liquid water potential temperature
    ! still missing the correction of order one
    !########################################################################
    subroutine Thermo_Anelastic_THETA_L(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        real(wp) Rdv_ov_Rd, Cdv_ov_Cd

        ! ###################################################################
        Rdv_ov_Rd = Rdv/Rd
        Cdv_ov_Cd = Cdv/Cd

#define locPi theta

        call Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, locPi)
        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)

        case (MIXT_TYPE_AIRVAPOR)
            theta = T*locPi**((1.0_wp + s(:, 2)*Rdv_ov_Rd)/(1.0_wp + s(:, 2)*Cdv_ov_Cd))

        case (MIXT_TYPE_AIRWATER)
            theta = T*locPi**((1.0_wp + s(:, 2)*Rdv_ov_Rd)/(1.0_wp + s(:, 2)*Cdv_ov_Cd))
            theta = theta*exp(-(Lv0 - T*Cvl)*s(:, 3)/T/(Cd + s(:, 2)*Cdv))

        end select

#undef locPi

        return
    end subroutine Thermo_Anelastic_THETA_L

    !########################################################################
    ! Equivalent potential temperature
    ! still missing the correction of order one
    !########################################################################
    subroutine Thermo_Anelastic_THETA_E(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        real(wp) Cdl_ov_Cd

        ! ###################################################################
        Cdl_ov_Cd = Cdl/Cd

#define locPi theta

        call Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, locPi)
        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)

        case (MIXT_TYPE_AIRVAPOR)
            theta = T*locPi**((1.0_wp - s(:, 2))/(1.0_wp + s(:, 2)*Cdl_ov_Cd))
            theta = theta*exp((Lv0 - T*Cvl)*s(:, 2)/T/(Cd + s(:, 2)*Cdl))

        case (MIXT_TYPE_AIRWATER)
            theta = T*locPi**((1.0_wp - s(:, 2))/(1.0_wp + s(:, 2)*Cdl_ov_Cd))
            theta = theta*exp((Lv0 - T*Cvl)*(s(:, 2) - s(:, 3))/T/(Cd + s(:, 2)*Cdl))

        end select

#undef locPi

        return
    end subroutine Thermo_Anelastic_THETA_E

    !########################################################################
    !########################################################################
    ! Frozen lapse rate; in unsaturated conditions, this is the unsaturated lapse rate
    subroutine Thermo_Anelastic_LAPSE_FR(nx, ny, nz, s, lapse)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: lapse(nx*ny*nz)

    ! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_AIR)
            lapse(:) = GRATIO*scaleheightinv

        case (MIXT_TYPE_AIRVAPOR)
            lapse(:) = GRATIO*scaleheightinv/(Cd + s(:, 2)*Cdv)

        case (MIXT_TYPE_AIRWATER)
            lapse(:) = GRATIO*scaleheightinv/(Cd + s(:, 2)*Cdv + s(:, 3)*Cvl)

        end select

        return
    end subroutine Thermo_Anelastic_LAPSE_FR

    !########################################################################
    !########################################################################
    ! Equilibrium lapse rate
    subroutine Thermo_Anelastic_LAPSE_EQU(nx, ny, nz, s, lapse, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: lapse(nx*ny*nz)
        real(wp), intent(inout) :: T(nx*ny*nz)

! -------------------------------------------------------------------
        real(wp) one_ov_Rd, one_ov_Rv, Rv_ov_Rd

! ###################################################################

        select case (imixture)
        case (MIXT_TYPE_AIR)
            call Thermo_Anelastic_LAPSE_FR(nx, ny, nz, s, lapse)

        case (MIXT_TYPE_AIRVAPOR)
            call Thermo_Anelastic_LAPSE_FR(nx, ny, nz, s, lapse)

        case (MIXT_TYPE_AIRWATER)

#define psat lapse
#define p_ov_psat lapse
#define qv_ov_qd lapse

            call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)

            call Thermo_Psat_Polynomial(nx*ny*nz, T, psat)

            p_ov_psat = 1.0_wp/psat
            call Thermo_Anelastic_WEIGHT_INPLACE(nx, ny, nz, pbackground, p_ov_psat)

            qv_ov_qd = rd_ov_rv/(p_ov_psat - 1.0_wp)

            one_ov_Rd = 1.0_wp/(Rd*GRATIO)
            one_ov_Rv = 1.0_wp/(Rv*GRATIO)
            Rv_ov_Rd = Rv/Rd
            lapse = (1.0_wp + qv_ov_qd*(Lv0 - T*Cvl)*one_ov_Rd/T) &
                    /(Cd + s(:, 2)*Cdl - qv_ov_qd*(1.0 - s(:, 2))*Cvl &
                      + qv_ov_qd*(1.0_wp - s(:, 2))*(1.0_wp + qv_ov_qd*Rv_ov_Rd)*(Lv0 - T*Cvl)**2.0_wp*one_ov_Rv/(T*T)) &
                    *GRATIO*scaleheightinv

#undef psat
#undef p_ov_psat
#undef qv_ov_qd

        end select

        return
    end subroutine Thermo_Anelastic_LAPSE_EQU

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_VAPOR_PRESSURE(nx, ny, nz, s, pv)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: pv(nx*ny*nz)

        ! ###################################################################
#define p pv
        call Thermo_Anelastic_PRESSURE(nx, ny, nz, p)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            pv = 0.0_wp

        case (MIXT_TYPE_AIRVAPOR)
            pv = s(:, 2)*Rv/(Rd + s(:, 2)*Rdv)*p

        case (MIXT_TYPE_AIRWATER)
            pv = (s(:, 2) - s(:, 3))*Rv/(Rd + s(:, 2)*Rdv - s(:, 3)*Rv)*p

        end select

#undef p

        return
    end subroutine Thermo_Anelastic_VAPOR_PRESSURE

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_RELATIVEHUMIDITY(nx, ny, nz, s, rh, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: rh(nx*ny*nz)
        real(wp), intent(inout) :: T(nx*ny*nz)

        ! ###################################################################
#define pv T
#define psat rh
        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)
        call Thermo_Psat_Polynomial(nx*ny*nz, T, psat)

        call Thermo_Anelastic_VAPOR_PRESSURE(nx, ny, nz, s, pv)

        rh = pv/psat*100.0_wp

#undef pv
#undef psat

        return
    end subroutine Thermo_Anelastic_RELATIVEHUMIDITY

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_DEWPOINT(nx, ny, nz, s, dpvdy, Td, Lapse)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *), dpvdy(nx*ny*nz)
        real(wp), intent(out) :: lapse(nx*ny*nz), Td(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) inr, nrmax
        real(wp) dpsat

        ! ###################################################################
        if (imixture == MIXT_TYPE_AIR) return

        ! maximum number of iterations in Newton-Raphson
        nrmax = 5

#define pv Lapse

        call Thermo_Anelastic_VAPOR_PRESSURE(nx, ny, nz, s, pv)

        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, Td)    ! Using actual temperature as initial condition

        ! If liquid, we should use qsat instead of ql (s(:,3)) because the smoothing function imposes
        ! an exponentially small value, but nonzero (?)

        do ij = 1, nx*ny*nz
            do inr = 1, nrmax                                   ! executing Newton-Raphson
                psat = THERMO_PSAT(NPSAT); dpsat = 0.0_wp
                do ipsat = NPSAT - 1, 1, -1
                    psat = psat*Td(ij) + THERMO_PSAT(ipsat)
                    dpsat = dpsat*Td(ij) + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                end do
                Td(ij) = Td(ij) - (psat - pv(ij))/dpsat
            end do

            Lapse(ij) = -dpvdy(ij)/dpsat

        end do

#undef pv

        return
    end subroutine Thermo_Anelastic_DEWPOINT

    !########################################################################
    !########################################################################
    ! Just to check what the effect of using a wrong cp would be
    subroutine Thermo_Anelastic_STATIC_CONSTANTCP(nx, ny, nz, s, result)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: result(nx*ny*nz)

    ! ###################################################################
#define T result
        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            result = s(:, 1)

        case (MIXT_TYPE_AIRVAPOR)
            result = s(:, 1) + (Cd - (Cd + s(:, 2)*Cdv))*T

        case (MIXT_TYPE_AIRWATER)
            result = s(:, 1) + (Cd - (Cd + s(:, 2)*Cdv + s(:, 3)*Cvl))*T

        end select

#undef T

        return
    end subroutine Thermo_Anelastic_STATIC_CONSTANTCP

! ! ###################################################################
! ! ###################################################################
!     subroutine Thermo_Anelastic_PH_RE(nx, ny, nz, s, e, p, wrk3d)
!         use THERMO_AIRWATER

!         integer(wi) nx, ny, nz
!         real(wp), intent(in) :: p(nx*ny*nz), e(nx*ny*nz)
!         real(wp), intent(inout) :: s(nx*ny*nz, *)
!         real(wp), intent(OUT) :: wrk3d(nx*ny*nz)

! ! -------------------------------------------------------------------
!         integer(wi) i, jk, iter, niter
!         real(wp) r_loc(1), en_loc(1), s_loc(2), dummy(1)
!         real(wp) p_loc(1), e_loc(1), t_loc(1)

! ! ###################################################################
!         niter = 5

!         s(:, 3) = 0.0_wp ! initialize, q_l=0

!         do iter = 1, niter ! iteration
! ! calculate density in wrk3d
!             call Thermo_Anelastic_DENSITY(nx, ny, nz, s, e, p, wrk3d)

! ! calculate energy
!             ij = 0
!             do jk = 0, ny*nz - 1
!                 is = mod(jk, ny) + 1
!                 P_LOC = pbackground(is)
!                 E_LOC = epbackground(is)

!                 do i = 1, nx
!                     ij = ij + 1

!                     r_loc = wrk3d(ij)
!                     s_loc(1) = s(ij, 2)
!                     s_loc(2) = s(ij, 3)
!                     en_loc = s(ij, 1) - E_LOC - CRATIO_INV*P_LOC/r_loc

! ! solve equilibrium (rho,e,q_i)
!                     call THERMO_AIRWATER_RE(1, s_loc, en_loc, r_loc, t_loc(:), dummy(:))

!                     s(ij, 3) = s_loc(2)

!                 end do
!             end do
!         end do

!         return
!     end subroutine Thermo_Anelastic_PH_RE

end module Thermo_Anelastic
