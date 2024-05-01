#include "dns_const.h"

module FI_RADIATION
    use TLAB_CONSTANTS, only: wp, wi, BCS_MAX, BCS_MIN
    use TLAB_VARS, only: imode_eqns, inb_scal_array
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_ARRAYS, only: wrk3d
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
    use OPR_ODES
    implicit none

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_LW_BULK1D_LIQUID = 1
    integer, parameter :: TYPE_LW_BULK1D = 2

    real(wp) :: sigma = 5.67037442e-8_wp ! W /m^2 /K

    public :: FI_RADIATION_INITIALIZE
    public :: FI_RADIATION_X

contains
!########################################################################
!########################################################################
    subroutine FI_RADIATION_INITIALIZE()
        ! in case nondimensional we need to adjust sigma

        return
    end subroutine FI_RADIATION_INITIALIZE

!########################################################################
!########################################################################
    subroutine LWBULK1_LIQUID(radiation, nx, ny, nz, g, a, source, flux)
        type(term_dt), intent(in) :: radiation
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: a(nx*nz, ny)                ! bulk absorption coefficent
        real(wp), intent(out) :: source(nx*nz, ny)
        real(wp), intent(out), optional :: flux(nx*nz, ny)

        target a, source, flux

! -----------------------------------------------------------------------
        integer(wi) j, nxy, nxz
        real(wp) f0, f1
        real(wp), pointer :: p_org(:, :) => null()
        real(wp), pointer :: p_tau(:, :) => null()
        real(wp), pointer :: p_source(:, :) => null()
        real(wp), pointer :: p_flux(:, :) => null()

! #######################################################################
        nxy = nx*ny     ! For transposition to make y direction the last one
        nxz = nx*nz

        if (nz == 1) then
            p_org => a
            p_tau(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            p_source => source
            if (present(flux)) p_flux => flux
        else
            p_org => source
            p_tau(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            if (present(flux)) then
                p_source => flux
                p_flux(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            else
                p_source(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            end if

#ifdef USE_ESSL
            call DGETMO(a, nxy, nxy, nz, p_org, nz)
#else
            call DNS_TRANSPOSE(a, nxy, nz, nxy, p_org, nz)
#endif
        end if

! ###################################################################
        ! Calculate (negative) optical thickness; integrating from the top, to be checked
        p_tau(:, ny) = 0.0_wp   ! boundary condition
        call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MAX)
        ! p_tau(:, 1) = 0.0_wp     ! boundary condition
        ! call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MIN)

! Calculate exp(-\tau(z, zmax)/\mu)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

! ###################################################################
! Calculate heating rate
        f0 = radiation%parameters(1)*radiation%parameters(2)
        f1 = radiation%parameters(3)*radiation%parameters(2)
        if (abs(radiation%parameters(3)) > 0.0_wp) then
            do j = ny, 1, -1
                p_source(:, j) = p_org(:, j)*(p_tau(:, j)*f0 &                       ! downward flux
                                              + p_tau(:, 1)/p_tau(:, j)*f1)       ! upward flux
            end do
        else
            p_source = p_org*p_tau*f0
        end if

        if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
            call DGETMO(p_source, nz, nz, nxy, source, nxy)
#else
            call DNS_TRANSPOSE(p_source, nz, nxy, nz, source, nxy)
#endif
        end if

! ###################################################################
! Calculate radiative flux, if necessary
        if (present(flux)) then
            f0 = -radiation%parameters(1)*radiation%parameters(2)
            f1 = radiation%parameters(3)*radiation%parameters(2)
            if (abs(radiation%parameters(3)) > 0.0_wp) then
                do j = ny, 1, -1
                    p_flux(:, j) = p_tau(:, j)*f0 &                       ! downward flux
                                   + p_tau(:, 1)/p_tau(:, j)*f1       ! upward flux
                end do
            else
                p_flux = p_tau*f0
            end if

            if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
                call DGETMO(p_flux, nz, nz, nxy, flux, nxy)
#else
                call DNS_TRANSPOSE(p_flux, nz, nxy, nz, flux, nxy)
#endif
            end if

        end if

! ###################################################################
        nullify (p_org, p_tau, p_source, p_flux)

        return
    end subroutine LWBULK1_LIQUID

!########################################################################
!########################################################################
    subroutine FI_RADIATION_X(radiation, nx, ny, nz, g, s, source, a, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: radiation
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)
        real(wp), intent(inout) :: a(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        ! -----------------------------------------------------------------------
        real(wp) delta_inv

        !########################################################################
        ! bulk absorption coefficient
        delta_inv = 1.0_wp/radiation%parameters(2)
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx, ny, nz, rbackground, s(:, radiation%scalar(1)), a)
            a = delta_inv*a
        else
            a = delta_inv*s(:, radiation%scalar(1))
        end if

        select case (radiation%type)
        case (TYPE_LW_BULK1D_LIQUID, EQNS_RAD_BULK1D_LOCAL)
            if (present(flux)) then
                call LWBULK1_LIQUID(radiation, nx, ny, nz, g, a, source, flux)
            else
                call LWBULK1_LIQUID(radiation, nx, ny, nz, g, a, source)
            end if

        case (TYPE_LW_BULK1D)
            ! to be done

        end select

    end subroutine FI_RADIATION_X

end module
