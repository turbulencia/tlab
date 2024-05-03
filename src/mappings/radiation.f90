#include "dns_const.h"
#include "dns_error.h"

module Radiation
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, BCS_MAX, BCS_MIN, efile, MAX_PROF
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_VARS, only: imode_eqns, inb_scal_array
    use TLAB_VARS, only: infrared
    use TLAB_ARRAYS, only: wrk3d
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    use THERMO_VARS, only: imixture
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
    use OPR_ODES
    implicit none

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_LW_BULK1D_LIQUID = 1
    integer, parameter :: TYPE_LW_BULK1D = 2

    real(wp) :: sigma = 5.67037442e-8_wp ! W /m^2 /K

    public :: Radiation_Initialize
    public :: Radiation_Infrared

contains
    ! ###################################################################
    ! ###################################################################
    subroutine RADIATION_READBLOCK(bakfile, inifile, block, var)
        character(len=*), intent(in) :: bakfile, inifile, block
        type(term_dt), intent(out) :: var

        character(len=512) sRes
        integer(wi) idummy

        ! -------------------------------------------------------------------
        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

        if (var%type == EQNS_NONE) then        ! backwards compatibility, to be removed
            call SCANINICHAR(bakfile, inifile, block, 'Type', 'None', sRes)
            if (trim(adjustl(sRes)) == 'none') then; var%type = TYPE_NONE
            else if (trim(adjustl(sRes)) == 'irbulk1dliquid') then; var%type = TYPE_LW_BULK1D_LIQUID
            else if (trim(adjustl(sRes)) == 'irbulk1d') then; var%type = TYPE_LW_BULK1D
            else if (trim(adjustl(sRes)) == 'bulk1dlocal') then; var%type = EQNS_RAD_BULK1D_LOCAL
            else if (trim(adjustl(sRes)) == 'bulk1dglobal') then; var%type = EQNS_RAD_BULK1D_GLOBAL
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Radiation option.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end if

        var%active = .false.
        if (var%type /= EQNS_NONE) then
            call SCANINIINT(bakfile, inifile, block, 'Scalar', '1', idummy)
            var%active(idummy) = .true.

            var%parameters(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, block, 'Parameters', '1.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, var%parameters)

        end if

        return
    end subroutine RADIATION_READBLOCK

!########################################################################
!########################################################################
    subroutine Radiation_Initialize()

        ! -------------------------------------------------------------------
        ! By default, transport and radiation are caused by last scalar
        infrared%scalar = inb_scal_array

        if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            if (infrared%type /= EQNS_NONE) then
                infrared%active(inb_scal_array) = .true. ! liquid
                infrared%active(inb_scal_array + 1) = .true. ! buoyancy
            end if

        end if

        ! backwards compatibility
        if (any([EQNS_RAD_BULK1D_LOCAL, EQNS_RAD_BULK1D_GLOBAL] == infrared%type)) then
            infrared%parameters(1) = infrared%parameters(1)*infrared%parameters(2)
            infrared%parameters(3) = infrared%parameters(3)*infrared%parameters(2)
            infrared%parameters(2) = 1.0_wp/infrared%parameters(2)
        end if

        ! -------------------------------------------------------------------
        ! in case nondimensional we need to adjust sigma

        return
    end subroutine Radiation_Initialize

!########################################################################
!########################################################################
    subroutine IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, a, source, flux)
        type(term_dt), intent(in) :: infrared
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
        f0 = infrared%parameters(1)
        f1 = infrared%parameters(3)
        if (abs(infrared%parameters(3)) > 0.0_wp) then
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
            f0 = -infrared%parameters(1)
            f1 = infrared%parameters(3)
            if (abs(infrared%parameters(3)) > 0.0_wp) then
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
    end subroutine IR_Bulk1D_Liquid

!########################################################################
!########################################################################
    subroutine Radiation_Infrared(infrared, nx, ny, nz, g, s, source, a, b, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)
        real(wp), intent(inout) :: a(nx*ny*nz)      ! space for bulk absorption coefficient
        real(wp), intent(inout) :: b(nx*ny*nz)      ! space for emission function
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        ! -----------------------------------------------------------------------
        real(wp) kappa, kappal, kappav

        !########################################################################
        select case (infrared%type)
        case (TYPE_LW_BULK1D_LIQUID, EQNS_RAD_BULK1D_LOCAL)
            ! bulk absorption coefficient
            kappa = infrared%parameters(2)
            a = kappa*s(:, infrared%scalar(1))
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, a)
            end if

            if (present(flux)) then
                call IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, a, source, flux)
            else
                call IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, a, source)
            end if

        case (TYPE_LW_BULK1D)
            ! bulk absorption coefficients and emission function
            kappal = infrared%parameters(2)
            kappav = infrared%parameters(3)
            a = kappal*s(:, infrared%scalar(1)) + kappav*(s(:,2) - s(:, infrared%scalar(1)))
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, a)

                call THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, s, wrk3d)
            else
                ! calculate temperature
            end if
            b = sigma*wrk3d**4./pi_wp

            ! if (present(flux)) then
            !     call IR_Bulk1D(infrared, nx, ny, nz, g, a, b, source, flux)
            ! else
            !     call IR_Bulk1D(infrared, nx, ny, nz, g, a, b, source)
            ! end if

        end select

    end subroutine Radiation_Infrared

end module
