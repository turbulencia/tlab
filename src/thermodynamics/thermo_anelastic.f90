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

module THERMO_ANELASTIC
    use TLab_Constants, only: wp, wi
    use Thermodynamics, only: imixture, GRATIO, scaleheight
    use Thermodynamics, only: THERMO_PSAT, NPSAT
    use Thermodynamics, only: Rv, Rd, Rdv, Cd, Cdv, Lv0, Ld, Ldv, Cvl, Cdl, Cl, rd_ov_rv, rd_ov_cd, PREF_1000
    implicit none
    private

    integer(wi) ij, i, j, jk, is, ipsat
    real(wp) RMEAN, P_LOC, E_LOC, T_LOC, R_LOC, R_LOC_INV, psat

    public :: THERMO_ANELASTIC_BUOYANCY
    public :: THERMO_ANELASTIC_DENSITY
    public :: THERMO_ANELASTIC_DEWPOINT
    public :: THERMO_ANELASTIC_LAPSE_EQU
    public :: THERMO_ANELASTIC_LAPSE_FR
    public :: THERMO_ANELASTIC_LWP
    public :: THERMO_ANELASTIC_QVEQU
    public :: THERMO_ANELASTIC_RELATIVEHUMIDITY
    public :: THERMO_ANELASTIC_STATIC_CONSTANTCP
    public :: THERMO_ANELASTIC_STATIC_L
    public :: THERMO_ANELASTIC_TEMPERATURE
    public :: THERMO_ANELASTIC_THETA
    public :: THERMO_ANELASTIC_THETA_E
    public :: THERMO_ANELASTIC_THETA_L
    public :: THERMO_ANELASTIC_THETA_V
    public :: THERMO_ANELASTIC_WEIGHT_ADD
    public :: THERMO_ANELASTIC_WEIGHT_INPLACE
    public :: THERMO_ANELASTIC_WEIGHT_OUTPLACE
    public :: THERMO_ANELASTIC_WEIGHT_SUBSTRACT
    public :: THERMO_ANELASTIC_PH               ! We could use THERMO_ANELASTIC_AIRWATER...
    ! public :: THERMO_ANELASTIC_PH_RE

    ! background, reference profiles
    real(wp), allocatable, public :: pbackground(:)                     ! Pressure background profile
    real(wp), allocatable, public :: tbackground(:)                     ! Temperature
    real(wp), allocatable, public :: rbackground(:), ribackground(:)    ! Density and its inverse
    real(wp), allocatable, public :: epbackground(:)                    ! Potential energy
    real(wp), allocatable, public :: sbackground(:, :)                  ! Scalars

contains

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, s, T)
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
    end subroutine THERMO_ANELASTIC_TEMPERATURE

!########################################################################
! Calculating h_l - h; very similar to the temperature routine
!########################################################################
    subroutine THERMO_ANELASTIC_STATIC_L(nx, ny, nz, s, result)
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
    end subroutine THERMO_ANELASTIC_STATIC_L

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_DENSITY(nx, ny, nz, s, rho)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: rho(nx*ny*nz)

! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    rho(ij) = P_LOC/T_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    RMEAN = Rd + s(ij, 2)*Rdv
                    rho(ij) = P_LOC/(RMEAN*T_LOC)
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    RMEAN = Rd + s(ij, 2)*Rdv - s(ij, 3)*Rv
                    rho(ij) = P_LOC/(RMEAN*T_LOC)
                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_DENSITY

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_BUOYANCY(nx, ny, nz, s, b)
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
    end subroutine THERMO_ANELASTIC_BUOYANCY

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_QVEQU(nx, ny, nz, s, T, qvequ)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: T(nx*ny*nz), qvequ(nx*ny*nz)

! ###################################################################
        if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T(ij) = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    psat = 0.0_wp
                    do ipsat = NPSAT, 1, -1
                        psat = psat*T(ij) + THERMO_PSAT(ipsat)
                    end do
                    qvequ(ij) = rd_ov_rv*(1.0_wp - s(ij, 2))/(P_LOC/psat - 1.0_wp)
                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_QVEQU

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_RELATIVEHUMIDITY(nx, ny, nz, s, T, rh)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: T(nx*ny*nz), rh(nx*ny*nz)

! ###################################################################
        if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1

                    T(ij) = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)

                    psat = 0.0_wp
                    do ipsat = NPSAT, 1, -1
                        psat = psat*T(ij) + THERMO_PSAT(ipsat)
                    end do
                    rh(ij) = (s(ij, 2) - s(ij, 3))*P_LOC/psat*Rv/(Rd + s(ij, 2)*Rdv - s(ij, 3)*Rv)
                    rh(ij) = rh(ij)*100.0_wp

                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_RELATIVEHUMIDITY

!########################################################################
! Dry potential temperature
!########################################################################
    subroutine THERMO_ANELASTIC_THETA(nx, ny, nz, s, theta)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz)

        real(wp) PI_LOC

        ! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                PI_LOC = (PREF_1000/pbackground(is))**rd_ov_cd
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    theta(ij) = T_LOC*PI_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                PI_LOC = (PREF_1000/pbackground(is))**rd_ov_cd
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    theta(ij) = T_LOC*PI_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                PI_LOC = (PREF_1000/pbackground(is))**rd_ov_cd
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    theta(ij) = T_LOC*PI_LOC
                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_THETA

!########################################################################
! Virtual Potential temperature
!########################################################################
    subroutine THERMO_ANELASTIC_THETA_V(nx, ny, nz, s, theta)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz)

        real(wp) PI_LOC

        ! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                PI_LOC = (PREF_1000/pbackground(is))**rd_ov_cd 
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    theta(ij) = T_LOC*PI_LOC 
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                PI_LOC = (PREF_1000/pbackground(is))**rd_ov_cd
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    theta(ij) = T_LOC*(Rd + s(ij, 2)*Rdv)/Rd*PI_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                PI_LOC = (PREF_1000/pbackground(is))**rd_ov_cd 
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    theta(ij) = T_LOC*(Rd + s(ij, 2)*Rdv - s(ij, 3)*Rv)/Rd*PI_LOC 
                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_THETA_V

!########################################################################
! Liquid water potential temperature
!########################################################################
    subroutine THERMO_ANELASTIC_THETA_L(nx, ny, nz, s, theta)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz)

        real(wp) kappa, Cp_loc, Lv

! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)/PREF_1000
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    theta(ij) = T_LOC/P_LOC**rd_ov_cd
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)/PREF_1000
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    kappa = (Rd + s(ij, 2)*Rdv)/(Cd + s(ij, 2)*Cdv)*GRATIO
                    theta(ij) = T_LOC/P_LOC**kappa
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)/PREF_1000
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1

                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    Cp_loc = 1.0_wp/(Cd + s(ij, 2)*Cdv)
                    kappa = (Rd + s(ij, 2)*Rdv)*Cp_loc*GRATIO
                    Lv = Lv0 - T_LOC*Cvl
                    theta(ij) = T_LOC/P_LOC**kappa*exp(-Lv*s(ij, 3)*Cp_loc/T_LOC)

                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_THETA_L

!########################################################################
! Equivalent potential temperature
!########################################################################
    subroutine THERMO_ANELASTIC_THETA_E(nx, ny, nz, s, theta)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz)

        real(wp) kappa, Cp_loc, Lv

! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)/PREF_1000
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    theta(ij) = T_LOC/P_LOC**rd_ov_cd
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)/PREF_1000
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1

                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv)
                    Cp_loc = 1.0_wp/(Cd + s(ij, 2)*Cdl)
                    kappa = Rd*(1.0_wp - s(ij, 2))*Cp_loc*GRATIO
                    Lv = Lv0 - T_LOC*Cvl
                    theta(ij) = T_LOC/P_LOC**kappa*exp(Lv*s(ij, 2)*Cp_loc/T_LOC)

                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)/PREF_1000
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1

                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    Cp_loc = 1.0_wp/(Cd + s(ij, 2)*Cdl)
                    kappa = Rd*(1.0_wp - s(ij, 2))*Cp_loc*GRATIO
                    Lv = Lv0 - T_LOC*Cvl
                    theta(ij) = T_LOC/P_LOC**kappa*exp(Lv*(s(ij, 2) - s(ij, 3))*Cp_loc/T_LOC)

                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_THETA_E

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_LAPSE_FR(nx, ny, nz, s, dTdy, lapse, frequency)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *), dTdy(nx*ny*nz)
        real(wp), intent(out) :: lapse(nx*ny*nz), frequency(nx*ny*nz)

! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            lapse(:) = GRATIO/scaleheight

            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    frequency(ij) = (lapse(ij) + dTdy(ij))/T_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            lapse(:) = GRATIO/scaleheight/(Cd + s(:, 2)*Cdv)

            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    frequency(ij) = (lapse(ij) + dTdy(ij))/T_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            lapse(:) = GRATIO/scaleheight/(Cd + s(:, 2)*Cdv + s(:, 3)*Cvl)

            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)
                    frequency(ij) = (lapse(ij) + dTdy(ij))/T_LOC
                end do

            end do

        end if

        return
    end subroutine THERMO_ANELASTIC_LAPSE_FR

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_LAPSE_EQU(nx, ny, nz, s, dTdy, dqldy, lapse, frequency)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *), dTdy(nx*ny*nz), dqldy(nx*ny*nz)
        real(wp), intent(out) :: lapse(nx*ny*nz), frequency(nx*ny*nz)

! -------------------------------------------------------------------
        real(wp) RT_INV, dpsat

        real(wp) Cp_loc, Lv, scaleheightinv, one_p_eps, qvequ, qsat, dummy

! ###################################################################
        scaleheightinv = GRATIO/scaleheight

        if (imixture == MIXT_TYPE_AIR) then
            lapse(:) = scaleheightinv

            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = s(ij, 1) - E_LOC
                    frequency(ij) = (lapse(ij) + dTdy(ij))/T_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            lapse(:) = scaleheightinv/(Cd + s(:, 2)*Cdv)

            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    frequency(ij) = (lapse(ij) + dTdy(ij))/T_LOC
                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)
                R_LOC = rbackground(is)

                RT_INV = R_LOC/(P_LOC)/scaleheight

                do i = 1, nx
                    ij = ij + 1

                    Cp_loc = Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl ! I need it below
                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/Cp_loc

                    psat = THERMO_PSAT(NPSAT); dpsat = 0.0_wp
                    do ipsat = NPSAT - 1, 1, -1
                        psat = psat*T_LOC + THERMO_PSAT(ipsat)
                        dpsat = dpsat*T_LOC + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                    end do
                    dummy = rd_ov_rv/(P_LOC/psat - 1.0_wp)
                    qsat = dummy/(1.0_wp + dummy)

! We cannot use ql directly (s(:,3)) because the smoothing function imposes
! an exponentially small value, but nonzero
                    if (qsat >= s(ij, 2)) then
                        lapse(ij) = scaleheightinv/(Cd + s(ij, 2)*Cdv)
                        frequency(ij) = (lapse(ij) + dTdy(ij))/T_LOC

                    else
                        qvequ = dummy*(1.0_wp - s(ij, 2))

                        one_p_eps = 1.0_wp/(1.0_wp - psat/P_LOC)
                        Lv = Lv0 - T_LOC*Cvl

                        lapse(ij) = scaleheightinv + qvequ*one_p_eps*Lv*RT_INV
                        lapse(ij) = lapse(ij)/(Cp_loc + qvequ*one_p_eps*Lv*dpsat/psat)

                        frequency(ij) = qvequ*(lapse(ij)*dpsat/psat - RT_INV)*one_p_eps - dqldy(ij)
                        frequency(ij) = (lapse(ij) + dTdy(ij))/T_LOC &
                                        - frequency(ij)*Rv/(Rd + s(ij, 2)*Rdv - s(ij, 3)*Rv)

                    end if

                end do

            end do

        end if

        return
    end subroutine THERMO_ANELASTIC_LAPSE_EQU

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_DEWPOINT(nx, ny, nz, s, Td, Lapse)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: lapse(nx*ny*nz), Td(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) inr, nrmax
        real(wp) dpsat, dummy, qsat, scaleheightinv

!  real(wp) NEWTONRAPHSON_ERROR

! ###################################################################
! maximum number of iterations in Newton-Raphson
        nrmax = 5

        scaleheightinv = 1.0_wp/scaleheight

        if (imixture == MIXT_TYPE_AIR) then

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)
                R_LOC = rbackground(is)

                do i = 1, nx
                    ij = ij + 1

! Using actual temperature as initial condition
                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)

! executing Newton-Raphson
                    do inr = 1, nrmax
                        psat = THERMO_PSAT(NPSAT); dpsat = 0.0_wp
                        do ipsat = NPSAT - 1, 1, -1
                            psat = psat*T_LOC + THERMO_PSAT(ipsat)
                            dpsat = dpsat*T_LOC + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                        end do
! we seek root of function:    psat -P_LOC *s(ij,2) *Rv /( Rd +s(ij,2) *Rdv )
                        T_LOC = T_LOC - (psat - P_LOC*s(ij, 2)*Rv/(Rd + s(ij, 2)*Rdv))/dpsat
                    end do
!           NEWTONRAPHSON_ERROR = MAX(NEWTONRAPHSON_ERROR,ABS(psat/dpsat)/T_LOC)
                    Td(ij) = T_LOC
                    Lapse(ij) = scaleheightinv*R_LOC/(P_LOC)*psat/dpsat

                end do
            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                P_LOC = pbackground(is)
                E_LOC = epbackground(is)
                R_LOC = rbackground(is)

                do i = 1, nx
                    ij = ij + 1

! Using actual temperature as initial condition
                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)

! We cannot use ql directly (s(:,3)) because the smoothing function imposes
! and exponentially small value, but nonzero
                    psat = THERMO_PSAT(NPSAT)
                    do ipsat = NPSAT - 1, 1, -1
                        psat = psat*T_LOC + THERMO_PSAT(ipsat)
                    end do
                    dummy = rd_ov_rv/(P_LOC/psat - 1.0_wp)
                    qsat = dummy/(1.0_wp + dummy)

                    if (qsat <= s(ij, 2)) then
                        Td(ij) = T_LOC

                    else
! executing Newton-Raphson
                        do inr = 1, nrmax
                            psat = THERMO_PSAT(NPSAT); dpsat = 0.0_wp
                            do ipsat = NPSAT - 1, 1, -1
                                psat = psat*T_LOC + THERMO_PSAT(ipsat)
                                dpsat = dpsat*T_LOC + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                            end do
!                                 psat -P_LOC *(s(ij,2)-s(ij,3)) *Rv /( Rd +s(ij,2) *Rdv -s(ij,3) *Rv )
! we seek root of function:       psat -P_LOC *s(ij,2) *Rv /( Rd +s(ij,2) *Rdv )
                            T_LOC = T_LOC - (psat - P_LOC*s(ij, 2)*Rv/(Rd + s(ij, 2)*Rdv))/dpsat
                        end do
                        ! NEWTONRAPHSON_ERROR = MAX(NEWTONRAPHSON_ERROR,ABS(psat/dpsat)/T_LOC)
                        ! print*,NEWTONRAPHSON_ERROR
                        Td(ij) = T_LOC
                        Lapse(ij) = scaleheightinv*R_LOC/(P_LOC)*psat/dpsat

                    end if

                end do
            end do

        end if

        return
    end subroutine THERMO_ANELASTIC_DEWPOINT

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, weight, a)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(*)
        real(wp), intent(inout) :: a(nx, ny*nz)

! ###################################################################
        do jk = 1, ny*nz
            j = mod(jk - 1, ny) + 1

            a(1:nx, jk) = a(1:nx, jk)*weight(j)

        end do

        return
    end subroutine THERMO_ANELASTIC_WEIGHT_INPLACE

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx, ny, nz, weight, a, b)
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
    end subroutine THERMO_ANELASTIC_WEIGHT_OUTPLACE

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_WEIGHT_ADD(nx, ny, nz, weight, a, b)
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
    end subroutine THERMO_ANELASTIC_WEIGHT_ADD

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_WEIGHT_SUBSTRACT(nx, ny, nz, weight, a, b)
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
    end subroutine THERMO_ANELASTIC_WEIGHT_SUBSTRACT

!########################################################################
!########################################################################
    subroutine THERMO_ANELASTIC_LWP(nx, ny, nz, g, ql, lwp, wrk1d, wrk3d)
        use TLab_Types, only: grid_dt
        use Integration, only: Int_Simpson

        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: ql(nx*nz, ny)
        real(wp), intent(out) :: lwp(nx, nz)
        real(wp), intent(INOUT) :: wrk1d(ny)
        real(wp), intent(INOUT) :: wrk3d(nx, ny, nz)

! -------------------------------------------------------------------
        integer(wi) k

! ###################################################################
        call THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx, ny, nz, rbackground, ql, wrk3d)

        do k = 1, nz
            do i = 1, nx
                do j = 1, ny
                    wrk1d(j) = wrk3d(i, j, k)
                end do
                lwp(i, k) = Int_Simpson(wrk1d(1:ny), g%nodes(1:ny))
            end do
        end do

        return
    end subroutine THERMO_ANELASTIC_LWP

!########################################################################
!########################################################################
! Just to check what the effect of using a wrong cp would be
    subroutine THERMO_ANELASTIC_STATIC_CONSTANTCP(nx, ny, nz, s, result)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: result(nx*ny*nz)

        real(wp) Lv

! ###################################################################
        if (imixture == MIXT_TYPE_AIR) then
            result = s(:, 1)

        else if (imixture == MIXT_TYPE_AIRVAPOR) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1

                    T_LOC = (s(ij, 1) - E_LOC)/(Cd + s(ij, 2)*Cdv)
                    result(ij) = T_LOC*Cd + E_LOC

                end do

            end do

        else if (imixture == MIXT_TYPE_AIRWATER) then
            ij = 0
            do jk = 0, ny*nz - 1
                is = mod(jk, ny) + 1
                E_LOC = epbackground(is)

                do i = 1, nx
                    ij = ij + 1

                    T_LOC = (s(ij, 1) - E_LOC + s(ij, 3)*Lv0)/(Cd + s(ij, 2)*Cdv + s(ij, 3)*Cvl)

                    Lv = Lv0 - T_LOC*Cvl
                    result(ij) = T_LOC*Cd + E_LOC - s(ij, 3)*Lv

                end do

            end do
        end if

        return
    end subroutine THERMO_ANELASTIC_STATIC_CONSTANTCP

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
    subroutine THERMO_ANELASTIC_PH(nx, ny, nz, s, h)
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
    end subroutine THERMO_ANELASTIC_PH

! ! ###################################################################
! ! ###################################################################
!     subroutine THERMO_ANELASTIC_PH_RE(nx, ny, nz, s, e, p, wrk3d)
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
!             call THERMO_ANELASTIC_DENSITY(nx, ny, nz, s, e, p, wrk3d)

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
!     end subroutine THERMO_ANELASTIC_PH_RE

end module THERMO_ANELASTIC
