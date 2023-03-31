#include "dns_const.h"
#include "dns_error.h"

! Implementation of thermal equation of state
! Default is NSP species represented by the first NSP-1 mass fractions in array s
! Airwater considers the special case where qt=qv+ql and ql instead of qv and ql are used as prognostic

! we use explicit-shape arrays in argements for the routines to be calleable with different array ranks

module THERMO_CALORIC
    use TLAB_CONSTANTS, only: wp, wi, MAX_NSP, efile
    use TLAB_PROCS
    use THERMO_VARS, only: imixture, THERMO_R, NSP
    use THERMO_VARS, only: gama0, CRATIO_INV, NCP, THERMO_AI, THERMO_TLIM
    use THERMO_VARS, only: Cd, Cdv, Cvl, Ld, Ldv, Lvl, Rd, Rdv, Rv
    implicit none
    private

    integer(wi) ij, is, icp, im
    real(wp) RMEAN, CPMEAN, CPMEAN_I, YMASS(MAX_NSP), ENTHALPY_I, dummy

    public :: THERMO_CALORIC_ENTHALPY
    public :: THERMO_CALORIC_ENERGY
    public :: THERMO_CALORIC_TEMPERATURE
    public :: THERMO_GAMMA
    public :: THERMO_CP

contains
    !########################################################################
    !# Calculate enthalpy from T and composition, h = \sum Y_ih_i
    !########################################################################
    subroutine THERMO_CALORIC_ENTHALPY(ijmax, s, T, h)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), T(ijmax)
        real(wp), intent(out) :: h(ijmax)

        ! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            h(:) = T(:)

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.
!     DO i = 1,nx*ny*nz
!#define MACRO_ZINPUT s(i,inb_scal)
!#include "dns_chem_mass.h"
!#define MACRO_TINPUT T(i)
!#include "dns_chem_enth.h"
!        h(i) = CH_H
!     ENDDO

        case (MIXT_TYPE_AIRWATER)   ! s(1,2) contains liquid mass fraction
            do ij = 1, ijmax
                h(ij) = (Cd + s(ij, 1)*Cdv + s(ij, 2)*Cvl)*T(ij) + (Ld + s(ij, 1)*Ldv + s(ij, 2)*Lvl)
            end do

        case default
            do ij = 1, ijmax
                ! pass species to YMASS vector
                YMASS(NSP) = 1.0_wp
                do is = 1, NSP - 1
                    YMASS(is) = s(ij, is)
                    YMASS(NSP) = YMASS(NSP) - YMASS(is)
                end do
                h(ij) = 0.0_wp
                do is = 1, NSP
                    if (T(ij) < THERMO_TLIM(3, is)) then
                        im = 2
                    else
                        im = 1
                    end if
                    ENTHALPY_I = 0.0_wp
                    do icp = NCP, 1, -1
                        ENTHALPY_I = ENTHALPY_I*T(ij) + THERMO_AI(icp, im, is)/real(icp, wp)
                    end do
                    ENTHALPY_I = ENTHALPY_I*T(ij) + THERMO_AI(6, im, is)
                    h(ij) = h(ij) + YMASS(is)*ENTHALPY_I
                end do

            end do

        end select

        return
    end subroutine THERMO_CALORIC_ENTHALPY

    !########################################################################
    !# Calculate energy from T and composition, e = \sum Y_ih_i - T\sum Y_i R_i
    !########################################################################
    subroutine THERMO_CALORIC_ENERGY(ijmax, s, T, e)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), T(ijmax)
        real(wp), intent(out) :: e(ijmax)

        ! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            dummy = 1.0/gama0 ! = 1.0_wp- RRATIO*CRATIO_INV
            e(:) = T(:)*dummy

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.
!     DO i = 1,nx*ny*nz
!#define MACRO_ZINPUT s(i,inb_scal)
!#include "dns_chem_mass.h"
!#define MACRO_TINPUT T(i)
!#include "dns_chem_enth.h"
!        e(i) = CH_H -...
!     ENDDO

        case (MIXT_TYPE_AIRWATER)   ! s(1,2) contains liquid mass fraction
            do ij = 1, ijmax
                e(ij) = (Cd + s(ij, 1)*Cdv + s(ij, 2)*Cvl)*T(ij) + (Ld + s(ij, 1)*Ldv + s(ij, 2)*Lvl)
                RMEAN = Rd + s(ij, 1)*Rdv - s(ij, 2)*Rv
                e(ij) = e(ij) - RMEAN*CRATIO_INV*T(ij)
            end do

        case default
            do ij = 1, ijmax
                ! pass species to YMASS vector
                YMASS(NSP) = 1.0_wp
                do is = 1, NSP - 1
                    YMASS(is) = s(ij, is)
                    YMASS(NSP) = YMASS(NSP) - YMASS(is)
                end do
                e(ij) = 0.0_wp
                RMEAN = 0.0_wp
                do is = 1, NSP
                    if (T(ij) < THERMO_TLIM(3, is)) then
                        im = 2
                    else
                        im = 1
                    end if
                    ENTHALPY_I = 0.0_wp
                    do icp = NCP, 1, -1
                        ENTHALPY_I = ENTHALPY_I*T(ij) + THERMO_AI(icp, im, is)/real(icp, wp)
                    end do
                    ENTHALPY_I = ENTHALPY_I*T(ij) + THERMO_AI(6, im, is)
                    e(ij) = e(ij) + YMASS(is)*ENTHALPY_I
                    RMEAN = RMEAN + YMASS(is)*THERMO_R(is)
                end do

                e(ij) = e(ij) - RMEAN*CRATIO_INV*T(ij)

            end do

        end select

        return
    end subroutine THERMO_CALORIC_ENERGY

    !########################################################################
    !# Calculate T from density, energy and composition, e=\sum Y_ih_i - T\sum Y_i R_i
    !# I need rho for the case in which the equilibrium composition is to be calculated.
    !########################################################################
    subroutine THERMO_CALORIC_TEMPERATURE(ijmax, s, e, rho, T, dqldqt)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), e(ijmax), rho(ijmax)
        real(wp), intent(out) :: T(ijmax)
        real(wp), intent(inout) :: dqldqt(ijmax)

        ! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            T(:) = gama0*e(:)

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.

        case (MIXT_TYPE_AIRWATER)   ! s(1,2) contains liquid mass fraction
            call THERMO_AIRWATER_RE(ijmax, s, e, rho, T, dqldqt)

        case default
            if (NCP == 1) then ! Cp linear with T. It assumes also only one temperature range
                do ij = 1, ijmax
                    CPMEAN = THERMO_AI(1, 1, NSP)               ! calculate heat capacity of mixture
                    T(ij) = e(ij) - THERMO_AI(6, 1, NSP)        ! substract formation energy of mixture
                    RMEAN = THERMO_R(NSP)                       ! calculate gas constant of mixture
                    do is = 1, NSP - 1
                        CPMEAN = CPMEAN + s(ij, is)*(THERMO_AI(1, 1, is) - THERMO_AI(1, 1, NSP))
                        T(ij) = T(ij) - s(ij, is)*(THERMO_AI(6, 1, is) - THERMO_AI(6, 1, NSP))
                        RMEAN = RMEAN + s(ij, is)*(THERMO_R(is) - THERMO_R(NSP))
                    end do
                    T(ij) = T(ij)/(CPMEAN - CRATIO_INV*RMEAN) ! solve for T; go from C_p to C_v
                end do

            else
                call TLAB_WRITE_ASCII(efile, 'THERMO_CALORIC_TEMPERATURE. General case undeveloped.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)

            end if

        end select

        return
    end subroutine THERMO_CALORIC_TEMPERATURE

    !########################################################################
    !# Calculate gama from T and composition.
    !########################################################################
    subroutine THERMO_GAMMA(ijmax, s, T, gama)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), T(ijmax)
        real(wp), intent(out) :: gama(ijmax)

        ! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            gama(:) = gama0

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.
            !     DO ij = 1,nx*ny*nz
            !#define MACRO_TINPUT T(ij)
            !#include "dns_chem_enth.h"
            !        gama(ij) = CPW/(CPW-GRATIO)
            !     ENDDO

        case (MIXT_TYPE_AIRWATER)   ! s(1,2) contains liquid mass fraction
            do ij = 1, ijmax
                CPMEAN = Cd + s(ij, 1)*Cdv + s(ij, 2)*Cvl
                RMEAN = Rd + s(ij, 1)*Rdv - s(ij, 2)*Rv
                gama(ij) = CPMEAN/(CPMEAN - RMEAN*CRATIO_INV)

            end do

        case default
            do ij = 1, ijmax
                ! pass species to YMASS vector
                YMASS(NSP) = 1.0_wp
                do is = 1, NSP - 1
                    YMASS(is) = s(ij, is)
                    YMASS(NSP) = YMASS(NSP) - YMASS(is)
                end do
                ! calculate cp and molecular weight of the mixture
                RMEAN = 0.0_wp
                CPMEAN = 0.0_wp
                do is = 1, NSP
                    if (T(ij) < THERMO_TLIM(3, is)) then
                        im = 2
                    else
                        im = 1
                    end if
                    CPMEAN_I = 0.0_wp
                    do icp = NCP, 1, -1
                        CPMEAN_I = CPMEAN_I*T(ij) + THERMO_AI(icp, im, is)
                    end do
                    CPMEAN = CPMEAN + YMASS(is)*CPMEAN_I
                    RMEAN = RMEAN + YMASS(is)*THERMO_R(is)
                end do

                gama(ij) = CPMEAN/(CPMEAN - RMEAN*CRATIO_INV)

            end do

        end select

        return
    end subroutine THERMO_GAMMA

    !########################################################################
    !# Calculate Cp from gamma and composition.
    !########################################################################
    subroutine THERMO_CP(ijmax, s, gama, cp)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: s(ijmax, *), gama(ijmax)
        real(wp), intent(out) :: cp(ijmax)

        ! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_NONE)
            cp(:) = 1.0_wp

        case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH, MIXT_TYPE_QUASIBS) ! Mass fractions defined by a conserved scalar Z.
!     DO ij = 1,ijmax
!#define MACRO_ZINPUT s(ij,inb_scal)
!#include "dns_chem_mass.h"
!        cp(ij) = gama(ij)*GRATIO/((gama(ij)-C_1_R)*WMEAN)
!     ENDDO

        case (MIXT_TYPE_AIRWATER)   ! s(1,2) contains liquid mass fraction
            do ij = 1, ijmax
                RMEAN = Rd + s(ij, 1)*Rdv - s(ij, 2)*Rv
                cp(ij) = gama(ij)*RMEAN*CRATIO_INV/(gama(ij) - 1.0_wp)
            end do

        case default
            do ij = 1, ijmax
                RMEAN = THERMO_R(NSP)
                do is = 1, NSP - 1
                    RMEAN = RMEAN + s(ij, is)*(THERMO_R(is) - THERMO_R(NSP))
                end do
                cp(ij) = gama(ij)*RMEAN*CRATIO_INV/(gama(ij) - 1.0_wp)
            end do

        end select

        return
    end subroutine THERMO_CP

end module THERMO_CALORIC
